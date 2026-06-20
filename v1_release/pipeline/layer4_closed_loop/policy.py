"""RL policy for the closed-loop selection layer.

Small MLP (~12k params) that scores each candidate given:
  - PCA(16)-reduced ESM3 embedding (16)
  - 4 engineered features (4)
  - World-model output: mean and variance for transduction / escape /
    inflammation (6)

Output: scalar `selection_score` per candidate.

Selection uses a Gumbel-top-k sampler (Plackett-Luce-style) over softmax-
normalized scores. This makes selection stochastic and gives a
differentiable log-prob proxy for REINFORCE updates.
"""

import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import config

POLICY_INPUT_DIM = 16 + 4 + 6  # PCA + engineered + GP (mean, var) × 3


class PolicyNetwork(nn.Module):
    """2-layer MLP. ~12k parameters."""

    def __init__(self, input_dim: int = POLICY_INPUT_DIM):
        super().__init__()
        self.input_dim = input_dim
        self.net = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.ReLU(),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """x: (N, input_dim) → scores (N,)."""
        return self.net(x).squeeze(-1)

    def num_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters())


def gumbel_top_k_sample(
    scores: torch.Tensor, k: int, tau: float, rng: np.random.Generator
) -> torch.Tensor:
    """Sample k distinct candidates without replacement via the
    Gumbel-top-k trick. Returns indices into `scores`.
    """
    n = scores.shape[0]
    if k >= n:
        return torch.arange(n)
    # Gumbel(0, 1) noise via inverse-CDF of uniform
    u = torch.tensor(rng.uniform(low=1e-9, high=1 - 1e-9, size=n), dtype=scores.dtype)
    g = -torch.log(-torch.log(u))
    perturbed = scores / tau + g
    return torch.topk(perturbed, k).indices


def selection_log_prob(scores: torch.Tensor, indices: torch.Tensor, tau: float) -> torch.Tensor:
    """Log-prob proxy of the chosen subset under softmax(scores/tau).

    Sum of per-pick log-softmax values — biased estimator (ignores
    without-replacement structure) but works in practice for REINFORCE.
    """
    log_p = torch.log_softmax(scores / tau, dim=0)
    return log_p[indices].sum()
