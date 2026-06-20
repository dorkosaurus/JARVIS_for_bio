"""Pretrain the RL policy on the wet-lab simulator.

Runs N campaigns of (10 cycles × 5 picks). Each campaign uses the same
80 candidate variants from `capsid_variants.fasta` but a fresh RNG, so the
simulator noise + Gumbel sampling are different. The policy learns a
*function* — score(features, GP output) — not specific candidate identities.

Output: `data/pretrained_policy.pt` (PyTorch state_dict + metadata).
"""

import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from sklearn.decomposition import PCA
from tqdm import trange

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config
from pipeline.layer3_world_model.world_model import (
    OUTPUTS,
    cosine,
    load_embedding,
)
from pipeline.layer4_closed_loop.closed_loop import (
    CandidatePool,
    run_campaign,
)
from pipeline.layer4_closed_loop.policy import PolicyNetwork


# ----------------------------- loaders ------------------------------------


def load_pool() -> CandidatePool:
    """Load the 80 generated capsid variants with embeddings + engineered features."""
    idx = pd.read_parquet(config.EMBEDDINGS_INDEX)
    aav7m8_hash = idx[idx.anchor_label == "AAV7m8"].sequence_hash.iloc[0]
    aav7m8_emb = load_embedding(aav7m8_hash)

    capsid_ids, embeddings, engineered, sim_inputs = [], [], [], []
    for rec in SeqIO.parse(config.CAPSID_VARIANTS_FASTA, "fasta"):
        cid = rec.id.split("|")[0]
        meta = dict(p.split("=") for p in rec.description.split("|")[1:] if "=" in p)
        has_ins = bool(int(meta["has_7mer"]))
        peptide = meta["insertion"] or None
        hamming = int(meta["hamming"])
        ins_len = len(peptide) if peptide else 0

        # Look up embedding
        row = idx[idx.capsid_id == cid].iloc[0]
        emb = load_embedding(row.sequence_hash)

        capsid_ids.append(cid)
        embeddings.append(emb)
        engineered.append(
            [int(has_ins), int(ins_len), int(hamming), cosine(emb, aav7m8_emb)]
        )
        sim_inputs.append(
            {
                "has_7mer_insertion": has_ins,
                "insertion_peptide": peptide,
                "insertion_length": ins_len,
                "hamming_to_aav2": hamming,
            }
        )

    return CandidatePool(
        capsid_ids=capsid_ids,
        embeddings=np.stack(embeddings),
        engineered=np.array(engineered, dtype=np.float64),
        sim_inputs=sim_inputs,
    )


def load_public_initial_obs() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Load the 45-row public dataset, deduped by capsid_id, with embeddings."""
    csvs = [
        "dalkara_2013.csv",
        "byrne_2020.csv",
        "kotterman_2015.csv",
        "reichel_bucher_inflammation.csv",
    ]
    frames = []
    for f in csvs:
        df = pd.read_csv(config.PUBLIC_DATA_DIR / f)
        for col in OUTPUTS:
            if col not in df.columns:
                df[col] = np.nan
        frames.append(df)
    df = pd.concat(frames, ignore_index=True)
    feat_cols = ["has_7mer_insertion", "insertion_length", "hamming_to_aav2"]
    df = df.groupby("capsid_id", as_index=False)[OUTPUTS + feat_cols].first()

    idx = pd.read_parquet(config.EMBEDDINGS_INDEX)
    by_cid = idx.set_index("capsid_id")["sequence_hash"].to_dict()
    aav7m8_hash = idx[idx.anchor_label == "AAV7m8"].sequence_hash.iloc[0]
    aav7m8_emb = load_embedding(aav7m8_hash)

    embeddings, engineered = [], []
    for _, row in df.iterrows():
        emb = load_embedding(by_cid[row.capsid_id])
        embeddings.append(emb)
        engineered.append(
            [
                int(row.has_7mer_insertion),
                int(row.insertion_length),
                int(row.hamming_to_aav2),
                cosine(emb, aav7m8_emb),
            ]
        )

    # Build initial observation rows (just the OUTPUTS + capsid_id)
    obs = df[["capsid_id"] + OUTPUTS].copy()
    obs["cycle"] = -1  # pretraining anchor
    obs["selection_strategy"] = "pretraining"
    obs["source"] = "public_literature"
    obs["meets_constraint"] = obs.inflammation_score < config.INFLAMMATION_THRESHOLD
    obs["is_simulated"] = False

    return obs, np.stack(embeddings), np.array(engineered, dtype=np.float64)


def load_pca_basis() -> PCA:
    """Fit PCA(16) on every unique embedding currently in the cache."""
    idx = pd.read_parquet(config.EMBEDDINGS_INDEX)
    hashes = idx.sequence_hash.unique()
    embeddings = np.stack([load_embedding(h) for h in hashes])
    pca = PCA(n_components=16)
    pca.fit(embeddings)
    return pca


# ----------------------------- pretraining --------------------------------


def main(
    n_campaigns: int = config.N_PRETRAIN_CAMPAIGNS,
    lr: float = 5e-3,
    tau: float = 0.5,
    seed: int = 0,
    baseline_ema_alpha: float = 0.1,
) -> None:
    print(f"Loading candidate pool, public data, PCA basis...")
    pool = load_pool()
    obs0, obs_emb, obs_eng = load_public_initial_obs()
    pca = load_pca_basis()
    print(f"  candidates:  {len(pool)}")
    print(f"  initial observations: {len(obs0)} rows")
    print(f"  PCA: 16 components, EVR={pca.explained_variance_ratio_.sum():.3f}")

    policy = PolicyNetwork()
    optimizer = torch.optim.Adam(policy.parameters(), lr=lr)
    print(f"\nPolicy: {policy.num_parameters():,} parameters")
    print(f"Pretraining: {n_campaigns} campaigns × {config.N_CYCLES} cycles × "
          f"{config.K_PER_CYCLE} picks  ({n_campaigns * config.N_CYCLES} GP refits)")

    rewards_per_campaign: list[float] = []
    per_cycle_reward_baseline = 0.0
    rng_master = np.random.default_rng(seed)
    t0 = time.time()
    for campaign_idx in trange(n_campaigns, desc="pretrain campaigns"):
        rng = np.random.default_rng(rng_master.integers(0, 2**31 - 1))
        result = run_campaign(
            pool=pool,
            initial_obs=obs0,
            initial_obs_embeddings=obs_emb,
            initial_obs_engineered=obs_eng,
            pca_basis=pca,
            policy=policy,
            optimizer=optimizer,
            rng=rng,
            n_cycles=config.N_CYCLES,
            k_per_cycle=config.K_PER_CYCLE,
            tau=tau,
            selection_strategy="rl_policy",
            fit_gp_hyperparameters=False,  # speed: skip MLL maximization
            reward_baseline=per_cycle_reward_baseline,
        )
        rewards = result["cycle_summary"]["reward"].to_numpy()
        total_reward = float(rewards.sum())
        rewards_per_campaign.append(total_reward)

        # Update EMA baseline using this campaign's mean per-cycle reward
        mean_per_cycle = float(rewards.mean())
        per_cycle_reward_baseline = (
            baseline_ema_alpha * mean_per_cycle
            + (1.0 - baseline_ema_alpha) * per_cycle_reward_baseline
        )

    t_elapsed = time.time() - t0
    print(f"\nDone in {t_elapsed:.1f}s ({t_elapsed/n_campaigns*1000:.1f} ms/campaign)")

    # Report reward trend
    window = max(1, n_campaigns // 10)
    early = float(np.mean(rewards_per_campaign[:window]))
    late = float(np.mean(rewards_per_campaign[-window:]))
    print(f"\nReward (HV improvement) trend:")
    print(f"  first {window} campaigns: mean total reward = {early:.4f}")
    print(f"  last  {window} campaigns: mean total reward = {late:.4f}")
    print(f"  improvement: {late - early:+.4f}")

    config.PRETRAINED_POLICY.parent.mkdir(parents=True, exist_ok=True)
    torch.save(
        {
            "state_dict": policy.state_dict(),
            "input_dim": policy.input_dim,
            "n_campaigns": n_campaigns,
            "lr": lr,
            "tau": tau,
            "rewards_per_campaign": rewards_per_campaign,
        },
        config.PRETRAINED_POLICY,
    )
    print(f"\nSaved policy to {config.PRETRAINED_POLICY}")


if __name__ == "__main__":
    main()
