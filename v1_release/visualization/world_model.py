"""World-model architecture diagram for the v1 AAV capsid optimization pipeline.

Matches the reference at repo root (world_model.png): 4 horizontal layers
with pastel backgrounds, colored boxes containing the components in each
layer, arrows showing data flow + closed-loop feedback.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config

# --- Layout constants -----------------------------------------------------

W, H = 14.0, 8.6
LAYER_H = 1.9       # height of each horizontal band
LAYER_GAP = 0.05    # gap between bands
LAYER_X0 = 0.55
LAYER_X1 = W - 0.25

# Each layer: (label, background color, box fill color, box text color)
LAYERS = [
    ("Layer 1 — Generative",        "#ede9fe", "#7c3aed", "white"),
    ("Layer 2 — Feature engineering","#dcfce7", "#15803d", "white"),
    ("Layer 3 — World model",       "#dbeafe", "#1d4ed8", "white"),
    ("Layer 4 — RL policy + wet lab","#fed7aa", "#c2410c", "white"),
]

# Box content per layer: list of (title, caption, optional marker)
BOXES = [
    [  # Layer 1
        ("Substitution variants",
         "1-6 AA swaps at\nVR-IV / V / VIII / IX",
         None),
        ("Insertion variants",
         "7-mer peptides @587\n(AAV.7m8 trick)",
         None),
        ("Stacked variants",
         "Insertion + 2-4 subs at\nnon-VIII VR loops",
         None),
    ],
    [  # Layer 2
        ("ESM3",
         "esmc-6b-2024-12 via Forge\n2560-d mean-pooled + pseudo-LL",
         None),
        ("Engineered features",
         "insertion flag, length,\nhamming, cosine-to-AAV.7m8",
         None),
        ("PCA(16)",
         "Dim reduction for the\nGP world model (95% EVR)",
         None),
    ],
    [  # Layer 3
        ("World model",
         "Multi-output GP (BoTorch)\nμ, σ² for transduction · escape · inflammation",
         None),
        ("Public data",
         "Dalkara · Byrne · Kotterman ·\nReichel (45 anchor rows)",
         None),
    ],
    [  # Layer 4
        ("RL policy",
         "MLP, 11.7k params. REINFORCE on\nPareto-hypervolume improvement.",
         None),
        ("Wet-lab simulator",
         "Load-bearing mock.\nCalibrated to literature.",
         "mock"),
    ],
]


def main(out_path: Path = config.OUTPUTS_DIR / "world_model.png") -> None:
    fig, ax = plt.subplots(figsize=(W, H))
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    # --- Title ---
    ax.text(W / 2, H - 0.25, "JARVIS-for-bio v1 — AAV capsid optimization (AMD)",
            ha="center", va="top", fontsize=15, fontweight="bold")
    ax.text(W / 2, H - 0.65,
            "Closed-loop world model: cached ESM3 embeddings + GP + RL policy + "
            "biologically grounded simulator",
            ha="center", va="top", fontsize=10, style="italic", color="#444")

    # --- Layers (top-to-bottom: Layer 1 on top) ---
    top_y = H - 1.1
    layer_centers = []
    for i, ((label, bg, fill, txt), boxes) in enumerate(zip(LAYERS, BOXES)):
        y1 = top_y - i * (LAYER_H + LAYER_GAP)
        y0 = y1 - LAYER_H
        layer_centers.append((y0, y1))

        # Background band
        ax.add_patch(Rectangle(
            (LAYER_X0, y0), LAYER_X1 - LAYER_X0, LAYER_H,
            facecolor=bg, edgecolor="none", zorder=1,
        ))
        # Layer label on the left (rotated)
        ax.text(LAYER_X0 - 0.18, (y0 + y1) / 2, label,
                ha="center", va="center", rotation=90,
                fontsize=10, fontweight="bold", color="#374151", zorder=2)

        # Boxes
        n = len(boxes)
        usable_w = LAYER_X1 - LAYER_X0 - 0.4
        gap = 0.30
        box_w = (usable_w - gap * (n - 1)) / n
        box_h = LAYER_H - 0.50
        for j, (title, caption, marker) in enumerate(boxes):
            x = LAYER_X0 + 0.20 + j * (box_w + gap)
            y = y0 + 0.25
            ax.add_patch(FancyBboxPatch(
                (x, y), box_w, box_h,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                facecolor=fill, edgecolor="black", linewidth=0.6,
                alpha=0.95, zorder=3,
            ))
            ax.text(x + box_w / 2, y + box_h - 0.18, title,
                    ha="center", va="top",
                    fontsize=11, fontweight="bold", color=txt, zorder=4)
            ax.text(x + box_w / 2, y + box_h - 0.55, caption,
                    ha="center", va="top",
                    fontsize=8.5, color=txt, zorder=4)
            if marker == "v2":
                ax.add_patch(FancyBboxPatch(
                    (x + box_w - 0.45, y + box_h - 0.30), 0.40, 0.25,
                    boxstyle="round,pad=0.01,rounding_size=0.04",
                    facecolor="white", edgecolor="black", linewidth=0.5, zorder=5,
                ))
                ax.text(x + box_w - 0.25, y + box_h - 0.175, "v2+",
                        ha="center", va="center",
                        fontsize=7, fontweight="bold", color="#1f2937", zorder=6)
            if marker == "mock":
                ax.add_patch(FancyBboxPatch(
                    (x + box_w - 0.55, y + box_h - 0.30), 0.50, 0.25,
                    boxstyle="round,pad=0.01,rounding_size=0.04",
                    facecolor="#fff7ed", edgecolor="#c2410c", linewidth=0.6, zorder=5,
                ))
                ax.text(x + box_w - 0.30, y + box_h - 0.175, "mock",
                        ha="center", va="center",
                        fontsize=7, fontweight="bold", color="#c2410c", zorder=6)

    # --- Down-arrows (data flow between layers) ---
    def down_arrow(x, y_top, y_bot, label=None, color="#374151"):
        a = FancyArrowPatch(
            (x, y_top), (x, y_bot),
            arrowstyle="-|>", mutation_scale=14,
            color=color, linewidth=1.4, zorder=7,
        )
        ax.add_patch(a)
        if label:
            ax.text(x + 0.12, (y_top + y_bot) / 2, label,
                    fontsize=8, color=color, va="center")

    # L1 -> L2
    down_arrow(W / 2, layer_centers[0][0] - 0.01, layer_centers[1][1] + 0.01,
               label="VP1 sequences")
    # L2 -> L3
    down_arrow(W / 2, layer_centers[1][0] - 0.01, layer_centers[2][1] + 0.01,
               label="features")
    # L3 -> L4
    down_arrow(W / 2, layer_centers[2][0] - 0.01, layer_centers[3][1] + 0.01,
               label="μ, σ² per output")

    # --- Feedback loop (right side, dashed green) ---
    fb_x = LAYER_X1 - 0.30
    a = FancyArrowPatch(
        (fb_x, layer_centers[3][1] + 0.05),
        (fb_x, layer_centers[2][0] - 0.02),
        arrowstyle="-|>", mutation_scale=14,
        color="#15803d", linewidth=1.6, linestyle=(0, (4, 3)), zorder=7,
        connectionstyle="arc3,rad=0",
    )
    ax.add_patch(a)
    ax.text(fb_x - 0.15, (layer_centers[3][1] + layer_centers[2][0]) / 2,
            "observations\nrefit GP + REINFORCE",
            fontsize=8, color="#15803d", ha="right", va="center", fontweight="bold")

    # --- Honesty footer ---
    fig.text(0.99, 0.01,
             "v1 release · pretrained policy + cached ESM3 embeddings are pre-computed indices · "
             "simulator is the load-bearing mock",
             ha="right", fontsize=7, color="#666", style="italic")

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
