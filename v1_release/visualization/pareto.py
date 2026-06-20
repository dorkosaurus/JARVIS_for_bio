"""Headline Pareto-frontier plot for AMD AAV capsid optimization.

Reads `outputs/pareto_data.parquet` and produces `outputs/pareto_frontier.png`.
The renderer is intentionally thin: if a column in the parquet is wrong,
the figure is wrong.

Layout (matches the conceptual reference at repo root):
  - X axis: rpe_transduction
  - Y axis: neut_escape
  - Grey dots: all evaluated candidates
  - Orange circles: RL-policy selections (color-graded by cycle)
  - Red diamond: AAV2 wildtype anchor
  - Green diamond: AAV.7m8 anchor
  - Blue dashed line: Pareto frontier (RL policy + public, constraint-passing)
  - Reference dashed lines: min thresholds for transduction + escape
  - Shaded regions: safe therapeutic window (green); inflammation violation (red)
  - Honesty footer: simulator label
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config

MIN_TRANSDUCTION = 0.10
MIN_ESCAPE = 0.50


def main(
    parquet_path: Path = config.OUTPUTS_DIR / "pareto_data.parquet",
    out_path: Path = config.OUTPUTS_DIR / "pareto_frontier.png",
) -> None:
    df = pd.read_parquet(parquet_path)

    fig, ax = plt.subplots(figsize=(10, 7))

    # --- Background: inflammation-violation hatched region (full plane) ---
    # Drawn first so other layers cover it.
    viol = df.dropna(subset=["rpe_transduction", "neut_escape"])
    viol = viol[viol["meets_constraint"] == False]
    ax.scatter(
        viol["rpe_transduction"], viol["neut_escape"],
        s=30, c="red", alpha=0.18, marker="x",
        label=f"Inflammation ≥ {config.INFLAMMATION_THRESHOLD} (filtered)",
    )

    # --- All evaluated candidates (grey backdrop) ---
    plot_all = df.dropna(subset=["rpe_transduction", "neut_escape"])
    ax.scatter(
        plot_all["rpe_transduction"], plot_all["neut_escape"],
        s=18, c="lightgrey", alpha=0.55, label="All evaluated", zorder=2,
    )

    # --- Random baseline picks (subtle) ---
    rand = plot_all[plot_all["selection_strategy"] == "random_baseline"]
    ax.scatter(
        rand["rpe_transduction"], rand["neut_escape"],
        s=28, c="#999999", marker="s", alpha=0.65,
        label="Random baseline picks", zorder=3,
    )

    # --- RL picks, colored by cycle ---
    rl = plot_all[plot_all["selection_strategy"] == "rl_policy"].copy()
    if len(rl) > 0:
        sc = ax.scatter(
            rl["rpe_transduction"], rl["neut_escape"],
            s=55, c=rl["cycle"], cmap="Oranges", edgecolors="black",
            linewidths=0.6, alpha=0.95, label="RL policy picks", zorder=5,
        )
        cb = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label("RL cycle", fontsize=9)

    # --- Anchors ---
    anchors = df[df["is_anchor"] == 1].drop_duplicates("anchor_label")
    for _, a in anchors.iterrows():
        if pd.isna(a["rpe_transduction"]):
            continue
        color = "red" if a["anchor_label"] == "AAV2" else "green"
        ax.scatter(
            a["rpe_transduction"], a["neut_escape"],
            s=240, marker="D", c=color, edgecolors="black",
            linewidths=1.2, zorder=8,
            label=f"{a['anchor_label']} anchor",
        )
        ax.annotate(
            a["anchor_label"],
            (a["rpe_transduction"], a["neut_escape"]),
            xytext=(8, 8), textcoords="offset points",
            fontsize=10, fontweight="bold", zorder=9,
        )

    # --- Pareto frontier (RL + public, constraint-passing) ---
    front = df[
        (df["is_on_pareto_frontier"] == 1)
        & (df["selection_strategy"].isin(["rl_policy", "pretraining"]))
    ].dropna(subset=["rpe_transduction", "neut_escape"])
    front = front.sort_values("rpe_transduction")
    if len(front) > 1:
        ax.plot(
            front["rpe_transduction"], front["neut_escape"],
            color="#1f5fbf", linewidth=1.6, linestyle="--",
            label="Pareto frontier", zorder=6,
        )

    # --- Reference threshold lines ---
    ax.axvline(MIN_TRANSDUCTION, color="grey", linestyle=":", alpha=0.6, linewidth=0.9)
    ax.axhline(MIN_ESCAPE, color="grey", linestyle=":", alpha=0.6, linewidth=0.9)

    # --- Target zone box ---
    ax.add_patch(plt.Rectangle(
        (0.35, 0.60), 0.65, 0.40,
        facecolor="lightgreen", alpha=0.12, edgecolor="green",
        linewidth=1.0, linestyle="--", zorder=1,
    ))
    ax.text(0.36, 0.95, "target zone", color="green", fontsize=9,
            fontweight="bold", alpha=0.7)

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("RPE transduction efficiency  (intravitreal route)", fontsize=11)
    ax.set_ylabel("Neutralizing-antibody escape  (fraction of panel)", fontsize=11)
    fig.suptitle("World model is exploring a Pareto frontier",
                 fontsize=13, fontweight="bold", y=0.99)
    ax.set_title(
        "RL policy navigates two engineering objectives under an "
        "inflammation-safety constraint",
        fontsize=10, style="italic", pad=8,
    )
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    ax.grid(alpha=0.25)

    # Honesty footer
    fig.text(
        0.99, 0.01,
        "Simulated data — wet-lab outcomes from biologically grounded simulator v1.0. See README.",
        ha="right", fontsize=7, color="#666666", style="italic",
    )

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
