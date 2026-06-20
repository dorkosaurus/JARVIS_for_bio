"""Headline convergence plot: Pareto hypervolume over cycles, RL vs random.

Reads `outputs/closed_loop_summary.csv` and produces
`outputs/rl_vs_random.png`.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config


def main(
    csv_path: Path = config.OUTPUTS_DIR / "closed_loop_summary.csv",
    out_path: Path = config.OUTPUTS_DIR / "rl_vs_random.png",
) -> None:
    df = pd.read_csv(csv_path)

    rl = df[df["selection_strategy"] == "rl_policy"].sort_values("cycle")
    rand = df[df["selection_strategy"] == "random_baseline"].sort_values("cycle")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        rl["cycle"], rl["pareto_hypervolume"],
        marker="o", linewidth=2.0, color="#d97706",
        label="RL policy (pretrained on simulator)",
    )
    ax.plot(
        rand["cycle"], rand["pareto_hypervolume"],
        marker="s", linewidth=1.8, color="#6b7280",
        label="Random baseline",
    )

    ax.set_xlabel("Closed-loop cycle", fontsize=11)
    ax.set_ylabel("Pareto hypervolume", fontsize=11)
    fig.suptitle("RL policy outperforms random selection",
                 fontsize=13, fontweight="bold", y=0.99)
    ax.set_title(
        f"k = {config.K_PER_CYCLE} picks per cycle · {config.N_CYCLES} cycles · "
        f"inflammation < {config.INFLAMMATION_THRESHOLD}",
        fontsize=9, style="italic", pad=8,
    )
    ax.legend(loc="lower right", fontsize=10)
    ax.grid(alpha=0.3)
    ax.set_xticks(range(int(df["cycle"].min()), int(df["cycle"].max()) + 1))

    fig.text(
        0.99, 0.01,
        "Simulated data — wet-lab outcomes from biologically grounded simulator v1.0.",
        ha="right", fontsize=7, color="#666666", style="italic",
    )

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
