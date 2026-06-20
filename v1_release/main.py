"""End-to-end pipeline runner for the AMD AAV capsid optimization v1 release.

Assumes Phases 1-4 are complete (sequences fetched, embeddings cached,
public data seeded, policy pretrained). This script:

  1. Loads the candidate pool and the pretrained policy
  2. Runs one RL campaign + one random campaign (same seed, same pool)
  3. Writes everything to SQLite (data/results.db)
  4. Exports outputs/variants.fasta, outputs/pareto_data.parquet,
     outputs/closed_loop_summary.csv

To regenerate the figures: `python visualization/pareto.py`.
"""

import gc
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config
from pipeline import db
from pipeline.layer4_closed_loop.closed_loop import run_campaign
from pipeline.layer4_closed_loop.policy import PolicyNetwork
from scripts.pretrain_policy import (
    load_pca_basis,
    load_pool,
    load_public_initial_obs,
)


# ----------------------------- helpers ------------------------------------


def _pareto_mask_2d(points: np.ndarray) -> np.ndarray:
    """Boolean mask of non-dominated rows (maximize both axes)."""
    n = len(points)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        p = points[i]
        dominated_by = (points >= p).all(axis=1) & (points > p).any(axis=1)
        if dominated_by.any():
            keep[i] = False
    return keep


def build_candidates_table(pool, idx_df: pd.DataFrame) -> pd.DataFrame:
    """Assemble the `candidates` table rows from the FASTA + embedding index."""
    cap_to_hash = idx_df.set_index("capsid_id")["sequence_hash"].to_dict()
    cap_to_anchor = idx_df.set_index("capsid_id")["anchor_label"].to_dict()

    rows = []
    for rec in SeqIO.parse(config.CAPSID_VARIANTS_FASTA, "fasta"):
        cid = rec.id.split("|")[0]
        meta = dict(p.split("=") for p in rec.description.split("|")[1:] if "=" in p)
        has_ins = bool(int(meta["has_7mer"]))
        peptide = meta["insertion"] or None
        hamming = int(meta["hamming"])
        seq = str(rec.seq)
        sha = cap_to_hash.get(cid)
        rows.append({
            "candidate_id": cid,
            "capsid_id": cid,
            "vp1_sequence": seq,
            "has_7mer_insertion": int(has_ins),
            "insertion_peptide": peptide,
            "hamming_to_aav2": hamming,
            "capsid_embedding_path": str(config.EMBEDDINGS_DIR / f"{sha}.h5") if sha else None,
            "is_anchor": 0,
            "anchor_label": None,
        })

    # Add anchors + public-data candidates
    for cid, sha in cap_to_hash.items():
        if cid in {r["candidate_id"] for r in rows}:
            continue
        # Pull sequence from public FASTA if present, else placeholder
        seq, has_ins, peptide, hamming = _lookup_public_meta(cid)
        anchor_label = cap_to_anchor.get(cid)
        anchor_label = None if pd.isna(anchor_label) else anchor_label
        rows.append({
            "candidate_id": cid,
            "capsid_id": cid,
            "vp1_sequence": seq,
            "has_7mer_insertion": int(has_ins),
            "insertion_peptide": peptide,
            "hamming_to_aav2": hamming,
            "capsid_embedding_path": str(config.EMBEDDINGS_DIR / f"{sha}.h5"),
            "is_anchor": int(anchor_label is not None),
            "anchor_label": anchor_label,
        })

    return pd.DataFrame(rows)


def _lookup_public_meta(cid: str) -> tuple[str, bool, str | None, int]:
    """Find sequence + variant metadata for a public-data candidate."""
    public_fasta = config.SEQUENCES_DIR / "public_capsids.fasta"
    if public_fasta.exists():
        for rec in SeqIO.parse(public_fasta, "fasta"):
            if rec.id.split("|")[0] == cid:
                meta = dict(p.split("=") for p in rec.description.split("|")[1:] if "=" in p)
                return (
                    str(rec.seq),
                    bool(int(meta.get("has_7mer", "0"))),
                    meta.get("insertion") or None,
                    int(meta.get("hamming", "0")),
                )
    return ("", False, None, 0)


def public_results_rows(obs0: pd.DataFrame) -> pd.DataFrame:
    """Convert the public-initial-observation DataFrame to result rows."""
    df = obs0.copy()
    df["cycle"] = -1
    df["selection_strategy"] = "pretraining"
    df["is_on_pareto_frontier"] = False  # set later
    if "source_version" not in df.columns:
        df["source_version"] = "literature"
    # Rows missing inflammation can't be evaluated against the constraint.
    df["meets_constraint"] = df["inflammation_score"].apply(
        lambda x: None if pd.isna(x) else bool(x < config.INFLAMMATION_THRESHOLD)
    )
    df["candidate_id"] = df["capsid_id"]
    return df[[
        "candidate_id", "cycle", "selection_strategy",
        "rpe_transduction", "neut_escape", "inflammation_score",
        "meets_constraint", "is_on_pareto_frontier",
        "source", "source_version", "is_simulated",
    ]]


def campaign_results_rows(observations: pd.DataFrame, strategy: str) -> pd.DataFrame:
    """Convert one campaign's observations into result rows for SQLite.

    Public-anchor rows already exist; filter to only new (sim-sourced) rows.
    """
    df = observations[observations["source"] == "wet_lab_simulator"].copy()
    df["candidate_id"] = df["capsid_id"]
    df["selection_strategy"] = strategy
    df["is_on_pareto_frontier"] = False  # set later
    return df[[
        "candidate_id", "cycle", "selection_strategy",
        "rpe_transduction", "neut_escape", "inflammation_score",
        "meets_constraint", "is_on_pareto_frontier",
        "source", "source_version", "is_simulated",
    ]]


def cycle_summary_rows(observations: pd.DataFrame, hv_history: list[float], strategy: str) -> pd.DataFrame:
    """Build the cycle_summary rows for one campaign."""
    sim = observations[observations["source"] == "wet_lab_simulator"]
    rows = []
    for c in sorted(sim["cycle"].unique()):
        cyc = sim[sim["cycle"] == c]
        n_violations = int((~cyc["meets_constraint"]).sum())
        points = cyc[["rpe_transduction", "neut_escape"]].to_numpy()
        rows.append({
            "cycle": int(c),
            "selection_strategy": strategy,
            "n_candidates_tested": len(cyc),
            "n_constraint_violations": n_violations,
            "best_rpe_transduction": float(cyc["rpe_transduction"].max()),
            "best_neut_escape": float(cyc["neut_escape"].max()),
            "pareto_frontier_size": int(_pareto_mask_2d(points).sum()) if len(points) else 0,
            "pareto_hypervolume": float(hv_history[int(c) + 1]),
        })
    return pd.DataFrame(rows)


def stamp_pareto_frontier(results: pd.DataFrame) -> pd.DataFrame:
    """Mark constraint-passing, non-dominated rows as Pareto frontier.

    Computed per selection_strategy across the union of (public + that
    strategy's sim observations). Public rows are always evaluated; sim
    rows compete within their strategy.
    """
    out = results.copy()
    out["is_on_pareto_frontier"] = False
    public_mask = out["selection_strategy"] == "pretraining"
    for strat in out["selection_strategy"].unique():
        if strat == "pretraining":
            continue
        rel = out[(out["selection_strategy"] == strat) | public_mask].copy()
        # meets_constraint True (1) only; NaN-inflammation public rows excluded.
        rel = rel[rel["meets_constraint"] == True].copy()
        # Also require both axes are observed (drop NaN x or y).
        rel = rel.dropna(subset=["rpe_transduction", "neut_escape"])
        if len(rel) == 0:
            continue
        pts = rel[["rpe_transduction", "neut_escape"]].to_numpy()
        keep = _pareto_mask_2d(pts)
        winning_ids = rel.index[keep]
        out.loc[winning_ids, "is_on_pareto_frontier"] = True
    return out


# ----------------------------- main ---------------------------------------


def main(seed: int = 2026) -> None:
    print("=== AMD AAV Capsid Optimization Pipeline (v1) ===\n")

    # --- Load everything ---
    print("[1/5] Loading pool + public data + PCA basis...")
    pool = load_pool()
    obs0, obs_emb, obs_eng = load_public_initial_obs()
    pca = load_pca_basis()
    idx_df = pd.read_parquet(config.EMBEDDINGS_INDEX)
    print(f"  candidate pool: {len(pool)}  public obs: {len(obs0)}")

    # --- Load pretrained policy ---
    print("[2/5] Loading pretrained policy...")
    ckpt = torch.load(config.PRETRAINED_POLICY, weights_only=False)
    policy = PolicyNetwork(input_dim=ckpt["input_dim"])
    policy.load_state_dict(ckpt["state_dict"])
    policy.eval()
    for p in policy.parameters():
        p.requires_grad_(False)
    print(f"  policy: {policy.num_parameters():,} params, trained {ckpt['n_campaigns']} campaigns")

    # --- Run campaigns ---
    print(f"[3/5] Running RL + random campaigns ({config.N_CYCLES} cycles, k={config.K_PER_CYCLE})...")
    t0 = time.time()
    with torch.no_grad():
        rl_result = run_campaign(
            pool=pool, initial_obs=obs0,
            initial_obs_embeddings=obs_emb, initial_obs_engineered=obs_eng,
            pca_basis=pca, policy=policy, optimizer=None,
            rng=np.random.default_rng(seed),
            n_cycles=config.N_CYCLES, k_per_cycle=config.K_PER_CYCLE,
            tau=ckpt["tau"], selection_strategy="rl_policy",
            fit_gp_hyperparameters=False, reward_baseline=0.0,
        )
        gc.collect()
        rand_result = run_campaign(
            pool=pool, initial_obs=obs0,
            initial_obs_embeddings=obs_emb, initial_obs_engineered=obs_eng,
            pca_basis=pca, policy=None, optimizer=None,
            rng=np.random.default_rng(seed),
            n_cycles=config.N_CYCLES, k_per_cycle=config.K_PER_CYCLE,
            tau=ckpt["tau"], selection_strategy="random_baseline",
            fit_gp_hyperparameters=False, reward_baseline=0.0,
        )
    print(f"  ran in {time.time()-t0:.1f}s")
    print(f"  RL    final HV = {rl_result['hv_history'][-1]:.4f}")
    print(f"  random final HV = {rand_result['hv_history'][-1]:.4f}")

    # --- Write to SQLite ---
    print("[4/5] Writing to SQLite (data/results.db)...")
    conn = db.init_db(config.RESULTS_DB, fresh=True)

    candidates = build_candidates_table(pool, idx_df)
    db.write_candidates(conn, candidates)

    public_rows = public_results_rows(obs0)
    rl_rows = campaign_results_rows(rl_result["observations"], "rl_policy")
    rand_rows = campaign_results_rows(rand_result["observations"], "random_baseline")
    all_results = pd.concat([public_rows, rl_rows, rand_rows], ignore_index=True)
    all_results = stamp_pareto_frontier(all_results)
    db.write_results(conn, all_results)

    rl_cs = cycle_summary_rows(rl_result["observations"], rl_result["hv_history"], "rl_policy")
    rand_cs = cycle_summary_rows(rand_result["observations"], rand_result["hv_history"], "random_baseline")
    db.write_cycle_summary(conn, pd.concat([rl_cs, rand_cs], ignore_index=True))

    print(f"  candidates: {len(candidates)}  results: {len(all_results)}  cycle_rows: {len(rl_cs)+len(rand_cs)}")

    # --- Export outputs ---
    print("[5/5] Exporting outputs/ ...")
    config.OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    variants_out = config.OUTPUTS_DIR / "variants.fasta"
    with open(variants_out, "w") as out:
        for fasta in [config.AAV2_FASTA, config.AAV7M8_FASTA, config.CAPSID_VARIANTS_FASTA]:
            out.write(fasta.read_text())
    n_synth = sum(1 for _ in SeqIO.parse(variants_out, "fasta"))
    pareto_df = db.export_pareto_parquet(conn, config.OUTPUTS_DIR / "pareto_data.parquet")
    summary_df = db.export_cycle_summary_csv(conn, config.OUTPUTS_DIR / "closed_loop_summary.csv")
    print(f"  outputs/variants.fasta             ({n_synth} synthesizable capsids)")
    print(f"  outputs/pareto_data.parquet        ({len(pareto_df)} rows)")
    print(f"  outputs/closed_loop_summary.csv    ({len(summary_df)} cycle rows)")

    conn.close()
    print("\nDone. Render figures with:")
    print("  python visualization/pareto.py")
    print("  python visualization/rl_vs_random.py")


if __name__ == "__main__":
    main()
