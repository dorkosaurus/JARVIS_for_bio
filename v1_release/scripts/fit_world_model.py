"""Fit the GP world model on public pre-training data and report R²/RMSE.

Joins the 4 public CSVs, dedupes by capsid_id (anchors appear in multiple
sources), looks up each capsid's ESM3 embedding from the cache, computes
the 4 engineered features, fits the three GP heads, and reports training-set
fit quality.

Phase 3 sanity check. Not the live closed loop; just a fit-quality report.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config
from pipeline.layer3_world_model.world_model import (
    OUTPUTS,
    AMDCapsidWorldModel,
    cosine,
    load_all_embeddings,
    load_embedding,
)

PUBLIC_CSVS = [
    "dalkara_2013.csv",
    "byrne_2020.csv",
    "kotterman_2015.csv",
    "reichel_bucher_inflammation.csv",
]


def load_public_data() -> pd.DataFrame:
    frames = []
    for fname in PUBLIC_CSVS:
        df = pd.read_csv(config.PUBLIC_DATA_DIR / fname)
        for col in OUTPUTS:
            if col not in df.columns:
                df[col] = np.nan
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def dedupe_by_capsid(df: pd.DataFrame) -> pd.DataFrame:
    """Anchors (AAV2, AAV7m8) appear in multiple CSVs; merge so each
    capsid_id has one row with the first non-NaN value per column."""
    feat_cols = [
        "has_7mer_insertion",
        "insertion_length",
        "hamming_to_aav2",
        "vp1_sequence",
    ]
    keep_cols = OUTPUTS + feat_cols
    return df.groupby("capsid_id", as_index=False)[keep_cols].first()


def lookup_embeddings(capsid_ids: list[str]) -> dict[str, np.ndarray]:
    idx = pd.read_parquet(config.EMBEDDINGS_INDEX)
    by_cid = idx.set_index("capsid_id")["sequence_hash"].to_dict()
    out: dict[str, np.ndarray] = {}
    for cid in capsid_ids:
        if cid not in by_cid:
            raise KeyError(
                f"capsid_id {cid} missing from embedding index. "
                "Re-run pipeline/layer2_features/esm3_embedder.py."
            )
        out[cid] = load_embedding(by_cid[cid])
    return out


def build_engineered(
    df: pd.DataFrame, embeddings: dict[str, np.ndarray]
) -> np.ndarray:
    aav7m8 = embeddings["AAV7m8"]
    rows = []
    for _, row in df.iterrows():
        cos = cosine(embeddings[row.capsid_id], aav7m8)
        rows.append(
            [
                int(row.has_7mer_insertion),
                int(row.insertion_length),
                int(row.hamming_to_aav2),
                cos,
            ]
        )
    return np.array(rows, dtype=np.float64)


def report(df: pd.DataFrame, preds: dict[str, tuple[np.ndarray, np.ndarray]]) -> None:
    print("\nTraining-set fit:")
    print(f"  {'output':<22} {'n':>4}  {'R²':>7}  {'RMSE':>7}  {'mean_var':>9}")
    for out in OUTPUTS:
        if out not in preds:
            print(f"  {out:<22}  (no model)")
            continue
        y = df[out].to_numpy()
        mask = ~np.isnan(y)
        mean, var = preds[out]
        if mask.sum() < 2:
            continue
        ss_res = ((y[mask] - mean[mask]) ** 2).sum()
        ss_tot = ((y[mask] - y[mask].mean()) ** 2).sum()
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        rmse = float(np.sqrt(ss_res / mask.sum()))
        mean_var = float(var[mask].mean())
        print(f"  {out:<22} {int(mask.sum()):>4}  {r2:>7.3f}  {rmse:>7.3f}  {mean_var:>9.4f}")


def main() -> None:
    df = load_public_data()
    print(f"loaded {len(df)} public-data rows across {len(PUBLIC_CSVS)} CSVs")

    df = dedupe_by_capsid(df)
    print(f"deduped to {len(df)} unique capsid_ids")

    emb_dict = lookup_embeddings(df.capsid_id.tolist())
    embeddings = np.stack([emb_dict[cid] for cid in df.capsid_id], axis=0)
    engineered = build_engineered(df, emb_dict)
    print(f"embeddings:  {embeddings.shape}  engineered: {engineered.shape}")

    all_emb = load_all_embeddings()
    print(f"PCA fit basis: {len(all_emb)} unique embeddings in cache")

    wm = AMDCapsidWorldModel(pca_dim=16)
    wm.fit_pca(all_emb)
    evr = float(wm.pca.explained_variance_ratio_.sum())
    print(f"PCA: kept {wm.pca_dim} components, cumulative explained variance = {evr:.3f}")

    Y = df[OUTPUTS]
    wm.fit(embeddings, engineered, Y)
    print(f"world model fit; heads trained: {list(wm.models.keys())}")

    preds = wm.predict(embeddings, engineered)
    report(df, preds)


if __name__ == "__main__":
    main()
