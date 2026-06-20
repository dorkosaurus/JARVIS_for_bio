"""Multi-output GP world model for AMD AAV capsid optimization.

Three independent SingleTaskGPs (BoTorch). Inputs per candidate are
PCA-reduced ESM3 embedding + 4 engineered features. Outputs:
  - rpe_transduction
  - neut_escape
  - inflammation_score

Each head fits only on rows whose label is non-NaN, so partial-label
sources (Kotterman 2015: escape only; Reichel 2017: inflammation only)
contribute where they can.
"""

import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.models.transforms import Normalize, Standardize
from gpytorch.mlls import ExactMarginalLogLikelihood
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import config

OUTPUTS = ["rpe_transduction", "neut_escape", "inflammation_score"]


def load_embedding(sequence_hash: str) -> np.ndarray:
    with h5py.File(config.EMBEDDINGS_DIR / f"{sequence_hash}.h5", "r") as f:
        return np.array(f["embedding"])


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def load_all_embeddings() -> np.ndarray:
    """Return one row per unique sequence_hash in the cache."""
    idx = pd.read_parquet(config.EMBEDDINGS_INDEX)
    hashes = idx.sequence_hash.unique()
    return np.stack([load_embedding(h) for h in hashes])


class AMDCapsidWorldModel:
    """GP world model: PCA(2560-d embedding) + 4 engineered features.

    Three independent SingleTaskGPs in a dict, keyed by output name.
    """

    def __init__(self, pca_dim: int = 16):
        self.pca_dim = pca_dim
        self.feature_dim = pca_dim + 4
        self.pca: PCA | None = None
        self.models: dict[str, SingleTaskGP] = {}

    # ---- PCA basis (fit once on the full embedding cache) ----

    def fit_pca(self, embeddings: np.ndarray) -> None:
        self.pca = PCA(n_components=self.pca_dim)
        self.pca.fit(embeddings)

    def featurize(
        self, embeddings: np.ndarray, engineered: np.ndarray
    ) -> np.ndarray:
        if self.pca is None:
            raise RuntimeError("PCA not fit; call fit_pca() first.")
        reduced = self.pca.transform(embeddings)
        return np.concatenate([reduced, engineered], axis=1)

    # ---- GP heads ----

    def fit(
        self,
        embeddings: np.ndarray,
        engineered: np.ndarray,
        Y: pd.DataFrame,
        fit_hyperparameters: bool = True,
    ) -> None:
        """Fit the three GP heads on the provided training rows.

        embeddings: (N, 2560) raw ESM3 mean-pooled
        engineered: (N, 4) [has_7mer, insertion_length, hamming_to_aav2, cosine_to_aav7m8]
        Y: DataFrame with OUTPUTS columns (NaN allowed)
        fit_hyperparameters: if False, skip MLL maximization and use default
            kernel hyperparameters. Useful during policy pretraining where
            we need fast GP refits and don't need optimal kernel tuning.
        """
        if self.pca is None:
            self.fit_pca(embeddings)
        X_full = self.featurize(embeddings, engineered)

        self.models = {}
        for out in OUTPUTS:
            y = Y[out].to_numpy(dtype=np.float64)
            mask = ~np.isnan(y)
            if mask.sum() < 3:
                continue
            X_t = torch.tensor(X_full[mask], dtype=torch.float64)
            y_t = torch.tensor(y[mask], dtype=torch.float64).unsqueeze(-1)
            model = SingleTaskGP(
                X_t,
                y_t,
                input_transform=Normalize(d=X_t.shape[1]),
                outcome_transform=Standardize(m=1),
            )
            if fit_hyperparameters:
                mll = ExactMarginalLogLikelihood(model.likelihood, model)
                fit_gpytorch_mll(mll)
            self.models[out] = model

    def predict(
        self, embeddings: np.ndarray, engineered: np.ndarray
    ) -> dict[str, tuple[np.ndarray, np.ndarray]]:
        """Return {output_name: (mean_array, var_array)} for each output."""
        X_full = self.featurize(embeddings, engineered)
        X_t = torch.tensor(X_full, dtype=torch.float64)
        out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for name, model in self.models.items():
            model.eval()
            with torch.no_grad():
                posterior = model.posterior(X_t)
                mean = posterior.mean.squeeze(-1).cpu().numpy()
                var = posterior.variance.squeeze(-1).cpu().numpy()
            out[name] = (mean, var)
        return out
