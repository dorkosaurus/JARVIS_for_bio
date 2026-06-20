"""Head-to-head: pretrained RL policy vs random baseline.

Loads the freshly pretrained policy, runs N trials of (RL campaign + random
campaign) sharing the same RNG seed per trial, and reports mean Pareto
hypervolume after 10 cycles.
"""

import gc
import sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config
from pipeline.layer4_closed_loop.closed_loop import run_campaign
from pipeline.layer4_closed_loop.policy import PolicyNetwork
from scripts.pretrain_policy import (
    load_pca_basis,
    load_pool,
    load_public_initial_obs,
)


def main(n_trials: int = 3) -> None:
    pool = load_pool()
    obs0, obs_emb, obs_eng = load_public_initial_obs()
    pca = load_pca_basis()

    ckpt = torch.load(config.PRETRAINED_POLICY, weights_only=False)
    policy = PolicyNetwork(input_dim=ckpt["input_dim"])
    policy.load_state_dict(ckpt["state_dict"])
    policy.eval()
    for p in policy.parameters():
        p.requires_grad_(False)

    rl_finals, random_finals, initial_hvs = [], [], []
    for trial in range(n_trials):
        rng_rl = np.random.default_rng(1000 + trial)
        rng_rand = np.random.default_rng(1000 + trial)  # match seed

        with torch.no_grad():
            rl_result = run_campaign(
                pool=pool,
                initial_obs=obs0,
                initial_obs_embeddings=obs_emb,
                initial_obs_engineered=obs_eng,
                pca_basis=pca,
                policy=policy,
                optimizer=None,
                rng=rng_rl,
                n_cycles=config.N_CYCLES,
                k_per_cycle=config.K_PER_CYCLE,
                tau=ckpt["tau"],
                selection_strategy="rl_policy",
                fit_gp_hyperparameters=False,
                reward_baseline=0.0,
            )
            rl_hvs = list(rl_result["hv_history"])
            del rl_result
            gc.collect()

            rand_result = run_campaign(
                pool=pool,
                initial_obs=obs0,
                initial_obs_embeddings=obs_emb,
                initial_obs_engineered=obs_eng,
                pca_basis=pca,
                policy=None,
                optimizer=None,
                rng=rng_rand,
                n_cycles=config.N_CYCLES,
                k_per_cycle=config.K_PER_CYCLE,
                tau=ckpt["tau"],
                selection_strategy="random_baseline",
                fit_gp_hyperparameters=False,
                reward_baseline=0.0,
            )
            rand_hvs = list(rand_result["hv_history"])
            del rand_result
            gc.collect()
        initial_hv = rl_hvs[0]
        initial_hvs.append(initial_hv)
        rl_finals.append(rl_hvs[-1])
        random_finals.append(rand_hvs[-1])
        print(
            f"trial {trial}: initial_hv={initial_hv:.4f}  "
            f"rl_final={rl_hvs[-1]:.4f}  random_final={rand_hvs[-1]:.4f}  "
            f"Δ={rl_hvs[-1]-rand_hvs[-1]:+.4f}",
            flush=True,
        )

    print()
    print(f"mean: initial={np.mean(initial_hvs):.4f}  "
          f"rl={np.mean(rl_finals):.4f}  "
          f"random={np.mean(random_finals):.4f}  "
          f"Δ={np.mean(rl_finals)-np.mean(random_finals):+.4f}")
    wins = sum(1 for r, n in zip(rl_finals, random_finals) if r > n)
    print(f"RL beats random in {wins}/{n_trials} trials")


if __name__ == "__main__":
    main()
