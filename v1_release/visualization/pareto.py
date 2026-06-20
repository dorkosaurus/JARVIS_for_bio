"""Headline Pareto-frontier plot for AMD AAV capsid optimization.

Reads `outputs/pareto_data.parquet` and produces `outputs/pareto_frontier.png`.

Layout (matches the conceptual reference at repo root):
  - X axis: rpe_transduction
  - Y axis: neut_escape
  - Blue dashed curve: **theoretical** Pareto frontier (engineering trade-off
    envelope; passes through the target zone in the safe therapeutic window)
  - Grey dots: all evaluated candidates (live BELOW the curve)
  - Orange circles: RL-policy selections (color-graded by cycle)
  - Red diamond: AAV2 wildtype anchor (baseline)
  - Green diamond: AAV.7m8 anchor (published intravitreal best)
  - Reference dashed lines: min thresholds for transduction + escape
  - Shaded regions: safe therapeutic window (green); inflammation violation (red x's)
  - Target zone: small box on the Pareto curve, inside the safe window
  - Honesty footer: simulator label
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config

MIN_TRANSDUCTION = 0.10
MIN_ESCAPE = 0.50

# Target-zone box: small region on the theoretical Pareto curve, inside the
# safe therapeutic window. Coordinates chosen so the curve passes through.
TARGET_ZONE = dict(x0=0.38, x1=0.58, y0=0.62, y1=0.78)


def theoretical_pareto_curve() -> tuple[np.ndarray, np.ndarray]:
    """The engineering trade-off envelope, monotone-decreasing.

    Control points calibrated so the curve:
      - peaks near (0, 0.95) — extreme escape-only variants exist in principle
      - bends through the target zone at ~(0.45, 0.70) — where the trade-off
        crosses the safe therapeutic window
      - decays to (0.95, 0.10) — fitness cliff past ~85% transduction
    """
    control_x = np.array([0.00, 0.08, 0.20, 0.32, 0.45, 0.58, 0.72, 0.85, 0.95])
    control_y = np.array([0.95, 0.88, 0.80, 0.74, 0.69, 0.60, 0.42, 0.22, 0.08])
    interp = PchipInterpolator(control_x, control_y)
    x = np.linspace(0.0, 0.95, 200)
    y = interp(x)
    return x, y


def main(
    parquet_path: Path = config.OUTPUTS_DIR / "pareto_data.parquet",
    out_path: Path = config.OUTPUTS_DIR / "pareto_frontier.png",
) -> None:
    df = pd.read_parquet(parquet_path)

    fig, ax = plt.subplots(figsize=(10, 8.4))

    # --- Safe therapeutic window (green shading) ---
    ax.add_patch(plt.Rectangle(
        (MIN_TRANSDUCTION, MIN_ESCAPE),
        1.02 - MIN_TRANSDUCTION, 1.02 - MIN_ESCAPE,
        facecolor="#bbf7d0", alpha=0.18, edgecolor="none", zorder=0,
    ))
    ax.text(
        0.96, 0.96, "safe therapeutic window",
        color="#15803d", fontsize=9, fontweight="bold",
        ha="right", va="top", style="italic", alpha=0.75,
        transform=ax.transData,
    )

    # --- Theoretical Pareto frontier (curved trade-off envelope) ---
    cx, cy = theoretical_pareto_curve()
    ax.plot(
        cx, cy, color="#1f5fbf", linewidth=2.2, linestyle="--",
        label="Pareto frontier (engineering trade-off)", zorder=4,
    )

    # --- Target zone (on the Pareto curve, inside the safe window) ---
    tz = TARGET_ZONE
    ax.add_patch(plt.Rectangle(
        (tz["x0"], tz["y0"]), tz["x1"] - tz["x0"], tz["y1"] - tz["y0"],
        facecolor="#fef3c7", alpha=0.55, edgecolor="#a16207",
        linewidth=1.5, linestyle="-", zorder=5,
    ))
    ax.text(
        (tz["x0"] + tz["x1"]) / 2, tz["y1"] + 0.025,
        "target zone\n(therapeutic efficacy\nwithin safety constraint)",
        color="#92400e", fontsize=9, fontweight="bold", ha="center", va="bottom",
        zorder=10,
    )

    # --- Constraint violations (red x's, drawn behind all candidates) ---
    viol = df.dropna(subset=["rpe_transduction", "neut_escape"])
    viol = viol[viol["meets_constraint"] == False]
    if len(viol) > 0:
        ax.scatter(
            viol["rpe_transduction"], viol["neut_escape"],
            s=30, c="red", alpha=0.20, marker="x",
            label=f"Inflammation ≥ {config.INFLAMMATION_THRESHOLD} (filtered)",
            zorder=2,
        )

    plot_all = df.dropna(subset=["rpe_transduction", "neut_escape"])

    # --- Grey backdrop: all evaluated candidates ---
    ax.scatter(
        plot_all["rpe_transduction"], plot_all["neut_escape"],
        s=16, c="lightgrey", alpha=0.6,
        label="All evaluated candidates", zorder=3,
    )

    # --- Random baseline picks ---
    rand = plot_all[plot_all["selection_strategy"] == "random_baseline"]
    ax.scatter(
        rand["rpe_transduction"], rand["neut_escape"],
        s=28, c="#9ca3af", marker="s", alpha=0.7,
        label="Random baseline picks", zorder=6,
    )

    # --- RL picks, colored by cycle ---
    rl = plot_all[plot_all["selection_strategy"] == "rl_policy"].copy()
    if len(rl) > 0:
        sc = ax.scatter(
            rl["rpe_transduction"], rl["neut_escape"],
            s=60, c=rl["cycle"], cmap="Oranges",
            edgecolors="black", linewidths=0.7, alpha=0.95,
            label="RL policy picks", zorder=7, vmin=0, vmax=config.N_CYCLES - 1,
        )
        cb = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label("RL cycle", fontsize=9)

    # --- Anchors ---
    anchors = df[df["is_anchor"] == 1].drop_duplicates("anchor_label")
    for _, a in anchors.iterrows():
        if pd.isna(a["rpe_transduction"]):
            continue
        color = "#dc2626" if a["anchor_label"] == "AAV2" else "#16a34a"
        ax.scatter(
            a["rpe_transduction"], a["neut_escape"],
            s=240, marker="D", c=color, edgecolors="black",
            linewidths=1.2, zorder=9,
            label=f"{a['anchor_label']} (anchor)",
        )
        # Position label offset varies by anchor to avoid the legend / axis area
        if a["anchor_label"] == "AAV2":
            xy_off = (14, 12)
        else:  # AAV7m8
            xy_off = (12, 8)
        ax.annotate(
            a["anchor_label"],
            (a["rpe_transduction"], a["neut_escape"]),
            xytext=xy_off, textcoords="offset points",
            fontsize=11, fontweight="bold", zorder=11,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                      edgecolor="none", alpha=0.85),
        )

    # --- Reference threshold lines (labels placed top-right of intersections) ---
    ax.axvline(MIN_TRANSDUCTION, color="grey", linestyle=":", alpha=0.6, linewidth=0.9)
    ax.axhline(MIN_ESCAPE, color="grey", linestyle=":", alpha=0.6, linewidth=0.9)
    # Move threshold labels to the upper-left corner so they don't collide with
    # the AAV2 anchor or the legend at the bottom.
    ax.text(MIN_TRANSDUCTION + 0.005, 1.005, f"min transduction = {MIN_TRANSDUCTION}",
            fontsize=8, color="#555555", style="italic", va="top")
    ax.text(0.005, MIN_ESCAPE - 0.005, f"min escape = {MIN_ESCAPE}",
            fontsize=8, color="#555555", style="italic", va="top")

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
    # Legend below the plot so it doesn't collide with anchors or the
    # threshold labels in any corner. Multi-column for compactness.
    ax.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.10),
        fontsize=9, framealpha=0.95, ncol=4, borderaxespad=0.0,
    )
    ax.grid(alpha=0.2)

    # Honesty footer
    fig.text(
        0.99, 0.01,
        "Simulated data — wet-lab outcomes from biologically grounded simulator v1.0. "
        "Pareto curve is the engineering envelope, not derived from evaluated points.",
        ha="right", fontsize=7, color="#666666", style="italic",
    )

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
