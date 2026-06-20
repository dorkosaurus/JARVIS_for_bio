"""Closed-loop campaign driver.

One campaign = N cycles × k candidates per cycle. Each cycle:
  1. Refit the GP world model on accumulated observations.
  2. Score every untested candidate with the RL policy (or random baseline).
  3. Filter out candidates whose predicted inflammation >= constraint.
  4. Sample top-k constraint-passing candidates (Gumbel-top-k).
  5. Run them through the evidence source (simulator in v1, wet lab in v2+).
  6. Update Pareto frontier + hypervolume; REINFORCE policy step.

Same driver runs pretraining (policy=trainable, persist=False) and the
live demo (policy=loaded weights, persist=True).
"""

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from botorch.utils.multi_objective.hypervolume import Hypervolume

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import config
from pipeline.layer3_world_model.world_model import (
    OUTPUTS,
    AMDCapsidWorldModel,
)
from pipeline.layer4_closed_loop.policy import (
    POLICY_INPUT_DIM,
    PolicyNetwork,
    gumbel_top_k_sample,
    selection_log_prob,
)
from pipeline.layer4_closed_loop.wet_lab_simulator import simulate_wet_lab

REF_POINT = torch.tensor([0.0, 0.0], dtype=torch.float64)


# ----------------------------- data classes -------------------------------


@dataclass
class CandidatePool:
    """All candidate metadata, keyed by capsid_id."""

    capsid_ids: list[str]
    embeddings: np.ndarray            # (N, 2560)
    engineered: np.ndarray            # (N, 4): has_7mer, ins_len, hamming, cosine_to_aav7m8
    sim_inputs: list[dict]            # per-candidate sim inputs (has_7mer, peptide, ins_len, hamming)

    def __len__(self) -> int:
        return len(self.capsid_ids)


@dataclass
class Campaign:
    """Stateful container for one closed-loop run."""

    pool: CandidatePool
    initial_obs: pd.DataFrame         # capsid_id + OUTPUTS columns
    initial_obs_embeddings: np.ndarray
    initial_obs_engineered: np.ndarray
    pca_basis: object                 # sklearn PCA (already fit)
    rng: np.random.Generator
    constraint_threshold: float = config.INFLAMMATION_THRESHOLD
    fit_gp_hyperparameters: bool = True  # set False during fast pretraining

    observations: pd.DataFrame = field(default_factory=pd.DataFrame)
    observation_embeddings: np.ndarray = field(default_factory=lambda: np.empty((0, 2560)))
    observation_engineered: np.ndarray = field(default_factory=lambda: np.empty((0, 4)))
    tested_ids: set = field(default_factory=set)
    hv_history: list = field(default_factory=list)

    def __post_init__(self):
        self.observations = self.initial_obs.copy()
        self.observation_embeddings = self.initial_obs_embeddings.copy()
        self.observation_engineered = self.initial_obs_engineered.copy()
        self.hv_history = [self._compute_hv(self.observations)]

    # ----------- helpers -----------

    def _compute_hv(self, obs: pd.DataFrame) -> float:
        passing = obs[obs.inflammation_score < self.constraint_threshold]
        if passing.empty:
            return 0.0
        pts = passing[["rpe_transduction", "neut_escape"]].dropna().to_numpy()
        if len(pts) == 0:
            return 0.0
        # Hypervolume over the Pareto frontier of (transduction, escape) with ref (0,0)
        pareto_mask = _pareto_frontier_mask(pts)
        pareto = pts[pareto_mask]
        if len(pareto) == 0:
            return 0.0
        hv = Hypervolume(ref_point=REF_POINT)
        return float(hv.compute(torch.tensor(pareto, dtype=torch.float64)))

    def _refit_world_model(self) -> AMDCapsidWorldModel:
        wm = AMDCapsidWorldModel(pca_dim=16)
        wm.pca = self.pca_basis
        Y = self.observations[OUTPUTS]
        wm.fit(
            self.observation_embeddings,
            self.observation_engineered,
            Y,
            fit_hyperparameters=self.fit_gp_hyperparameters,
        )
        return wm

    def _untested_indices(self) -> np.ndarray:
        return np.array(
            [i for i, cid in enumerate(self.pool.capsid_ids) if cid not in self.tested_ids],
            dtype=int,
        )

    def _build_policy_input(
        self,
        untested_idx: np.ndarray,
        wm: AMDCapsidWorldModel,
    ) -> tuple[torch.Tensor, dict[str, tuple[np.ndarray, np.ndarray]]]:
        """Return policy input (N, 26) and dict of GP predictions per output."""
        emb = self.pool.embeddings[untested_idx]
        eng = self.pool.engineered[untested_idx]
        wm_x = wm.featurize(emb, eng)  # (N, 20)

        preds = wm.predict(emb, eng)  # {out: (mean(N,), var(N,))}
        means = np.stack([preds[o][0] for o in OUTPUTS], axis=1)  # (N, 3)
        vars_ = np.stack([preds[o][1] for o in OUTPUTS], axis=1)  # (N, 3)

        policy_x = np.concatenate([wm_x, means, vars_], axis=1)  # (N, 26)
        return torch.tensor(policy_x, dtype=torch.float32), preds

    # ----------- main step -----------

    def cycle(
        self,
        cycle_idx: int,
        policy: PolicyNetwork | None,
        k: int,
        tau: float,
        selection_strategy: str,
    ) -> dict:
        """Run one cycle. policy=None means random selection.

        Returns a dict with per-cycle info: picks, outcomes, reward, log_p,
        constraint_violations, hypervolume.
        """
        wm = self._refit_world_model()

        untested_idx = self._untested_indices()
        if len(untested_idx) < k:
            raise RuntimeError(
                f"Not enough untested candidates ({len(untested_idx)} < {k}) for cycle {cycle_idx}."
            )

        policy_x, preds = self._build_policy_input(untested_idx, wm)
        infl_means = preds["inflammation_score"][0]
        passing_mask = infl_means < self.constraint_threshold
        passing_local_idx = np.where(passing_mask)[0]

        if len(passing_local_idx) < k:
            # Constraint too tight — relax and pick lowest-inflammation candidates.
            order = np.argsort(infl_means)
            passing_local_idx = order[: max(k, 1)]

        if policy is None:
            # Random baseline
            scores_passing = torch.tensor(
                self.rng.standard_normal(len(passing_local_idx)), dtype=torch.float32
            )
            log_p = torch.tensor(0.0)
        else:
            scores_all = policy(policy_x)             # (N_untested,)
            scores_passing = scores_all[passing_local_idx]

        # Sample k by Gumbel-top-k over the passing subset
        if policy is None:
            top = torch.topk(scores_passing, k).indices.tolist()
            log_p = None
        else:
            picked_local = gumbel_top_k_sample(scores_passing, k, tau, self.rng)
            log_p = selection_log_prob(scores_passing, picked_local, tau)
            top = picked_local.tolist()

        # Map back to pool indices
        picked_passing_local = [int(passing_local_idx[i]) for i in top]
        picked_pool_idx = [int(untested_idx[i]) for i in picked_passing_local]
        picked_ids = [self.pool.capsid_ids[i] for i in picked_pool_idx]

        # Simulate
        outcomes = []
        n_violations = 0
        prev_hv = self.hv_history[-1]
        for pool_i in picked_pool_idx:
            sim_input = self.pool.sim_inputs[pool_i]
            out = simulate_wet_lab(
                has_7mer_insertion=sim_input["has_7mer_insertion"],
                insertion_peptide=sim_input["insertion_peptide"],
                insertion_length=sim_input["insertion_length"],
                hamming_to_aav2=sim_input["hamming_to_aav2"],
                rng=self.rng,
            )
            if out["inflammation_score"] >= self.constraint_threshold:
                n_violations += 1
            outcomes.append(out)

        # Append observations
        new_rows = []
        for cid, pool_i, out in zip(picked_ids, picked_pool_idx, outcomes):
            new_rows.append(
                {
                    "capsid_id": cid,
                    "cycle": cycle_idx,
                    "selection_strategy": selection_strategy,
                    "rpe_transduction": out["rpe_transduction"],
                    "neut_escape": out["neut_escape"],
                    "inflammation_score": out["inflammation_score"],
                    "meets_constraint": out["inflammation_score"] < self.constraint_threshold,
                    "source": "wet_lab_simulator",
                    "source_version": config.SIMULATOR_VERSION,
                    "is_simulated": True,
                }
            )
        new_obs_df = pd.DataFrame(new_rows)
        self.observations = pd.concat([self.observations, new_obs_df], ignore_index=True)
        self.observation_embeddings = np.concatenate(
            [self.observation_embeddings, self.pool.embeddings[picked_pool_idx]], axis=0
        )
        self.observation_engineered = np.concatenate(
            [self.observation_engineered, self.pool.engineered[picked_pool_idx]], axis=0
        )
        for cid in picked_ids:
            self.tested_ids.add(cid)

        # Hypervolume after picks
        new_hv = self._compute_hv(self.observations)
        self.hv_history.append(new_hv)
        reward = new_hv - prev_hv

        return {
            "cycle": cycle_idx,
            "picked_ids": picked_ids,
            "outcomes": outcomes,
            "reward": float(reward),
            "log_p": log_p,
            "prev_hv": float(prev_hv),
            "new_hv": float(new_hv),
            "n_violations": n_violations,
            "n_constraint_passing": int(passing_mask.sum()),
        }


# ----------------------------- helpers ------------------------------------


def _pareto_frontier_mask(points: np.ndarray) -> np.ndarray:
    """Boolean mask of non-dominated rows (maximize both columns).

    Row i is dominated by row j if points[j] >= points[i] (all dims) and
    > points[i] in at least one dim — i.e. j is at least as good
    everywhere and strictly better somewhere. Such i's are removed.
    """
    n = len(points)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        p = points[i]
        is_dominated_by = (
            (points >= p).all(axis=1) & (points > p).any(axis=1)
        )
        if is_dominated_by.any():
            keep[i] = False
    return keep


def run_campaign(
    pool: CandidatePool,
    initial_obs: pd.DataFrame,
    initial_obs_embeddings: np.ndarray,
    initial_obs_engineered: np.ndarray,
    pca_basis: object,
    policy: PolicyNetwork | None,
    optimizer: torch.optim.Optimizer | None,
    rng: np.random.Generator,
    n_cycles: int = config.N_CYCLES,
    k_per_cycle: int = config.K_PER_CYCLE,
    tau: float = 1.0,
    selection_strategy: str = "rl_policy",
    fit_gp_hyperparameters: bool = True,
    reward_baseline: float = 0.0,
) -> dict:
    """Run one campaign end-to-end. Returns observations + cycle summary.

    `reward_baseline` is subtracted from per-cycle reward before the REINFORCE
    gradient — reduces variance dramatically when caller maintains an EMA of
    past rewards.
    """
    campaign = Campaign(
        pool=pool,
        initial_obs=initial_obs,
        initial_obs_embeddings=initial_obs_embeddings,
        initial_obs_engineered=initial_obs_engineered,
        pca_basis=pca_basis,
        rng=rng,
        fit_gp_hyperparameters=fit_gp_hyperparameters,
    )

    cycle_summaries = []
    for c in range(n_cycles):
        info = campaign.cycle(
            cycle_idx=c,
            policy=policy,
            k=k_per_cycle,
            tau=tau,
            selection_strategy=selection_strategy,
        )

        # REINFORCE step with baseline-subtracted advantage
        if policy is not None and optimizer is not None and info["log_p"] is not None:
            advantage = info["reward"] - reward_baseline
            loss = -torch.tensor(advantage, dtype=torch.float32) * info["log_p"]
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        cycle_summaries.append(
            {
                "cycle": c,
                "selection_strategy": selection_strategy,
                "n_candidates_tested": k_per_cycle,
                "n_constraint_violations": info["n_violations"],
                "n_constraint_passing_predicted": info["n_constraint_passing"],
                "reward": info["reward"],
                "prev_hv": info["prev_hv"],
                "new_hv": info["new_hv"],
            }
        )

    summary_df = pd.DataFrame(cycle_summaries)
    return {
        "observations": campaign.observations.copy(),
        "cycle_summary": summary_df,
        "hv_history": list(campaign.hv_history),
        "tested_ids": set(campaign.tested_ids),
    }
