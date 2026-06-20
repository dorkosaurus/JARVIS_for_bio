"""World-model architecture diagram for the v1 AAV capsid optimization pipeline.

4 horizontal pastel layers (generative / features / world model / RL+wet lab).
Layer labels sit in a dedicated left margin (no collision with bands).
Down arrows live in dedicated inter-layer gaps with labels on white space.
Right-side feedback loop shows observations refitting the GP + REINFORCE.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config

# --- Layout constants -----------------------------------------------------

W, H = 16.0, 11.0
LAYER_H = 1.95           # height of each pastel band
LAYER_GAP = 0.65         # whitespace gap between bands (for arrows + labels)
LEFT_MARGIN = 1.70       # dedicated column for wrapped layer names
RIGHT_MARGIN = 1.80      # dedicated column for the feedback-loop arrow
LAYER_X0 = LEFT_MARGIN
LAYER_X1 = W - RIGHT_MARGIN

LAYERS = [
    ("Layer 1\nGenerative",          "#ede9fe", "#7c3aed", "white"),
    ("Layer 2\nFeature\nengineering","#dcfce7", "#15803d", "white"),
    ("Layer 3\nWorld\nmodel",        "#dbeafe", "#1d4ed8", "white"),
    ("Layer 4\nRL policy\n+ wet lab","#fed7aa", "#c2410c", "white"),
]

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
         "MLP, 11.7k params. REINFORCE\non Pareto-hypervolume improvement.",
         None),
        ("Wet-lab simulator",
         "Load-bearing mock.\nCalibrated to literature.",
         "mock"),
    ],
]

# Labels for the inter-layer down arrows
ARROW_LABELS = [
    "VP1 sequences",
    "ESM3 features",
    "μ, σ² per output",
]


def main(out_path: Path = config.OUTPUTS_DIR / "world_model.png") -> None:
    fig, ax = plt.subplots(figsize=(W, H))
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    # --- Title block --------------------------------------------------------
    ax.text(W / 2, H - 0.30,
            "JARVIS-for-bio v1 — AAV capsid optimization (AMD)",
            ha="center", va="top", fontsize=16, fontweight="bold")
    ax.text(W / 2, H - 0.80,
            "Closed-loop world model: cached ESM3 embeddings + GP + RL policy + "
            "biologically grounded simulator",
            ha="center", va="top", fontsize=10.5, style="italic", color="#444")

    # --- Layers (top-to-bottom) --------------------------------------------
    top_y = H - 1.35
    layer_y = []  # list of (y_bottom, y_top) for each layer band
    for i, ((label, bg, fill, txt), boxes) in enumerate(zip(LAYERS, BOXES)):
        y1 = top_y - i * (LAYER_H + LAYER_GAP)
        y0 = y1 - LAYER_H
        layer_y.append((y0, y1))

        # Background pastel band (only spans the inner column)
        ax.add_patch(Rectangle(
            (LAYER_X0, y0), LAYER_X1 - LAYER_X0, LAYER_H,
            facecolor=bg, edgecolor="none", zorder=1,
        ))

        # Horizontal multi-line layer label in the dedicated left margin
        # (no rotation — easier to read, fits naturally inside the band height).
        ax.text(LEFT_MARGIN / 2, (y0 + y1) / 2, label,
                ha="center", va="center",
                fontsize=11, fontweight="bold", color="#1f2937",
                linespacing=1.25, zorder=2)

        # Boxes inside the band
        n = len(boxes)
        usable_w = LAYER_X1 - LAYER_X0 - 0.5
        inner_gap = 0.32
        box_w = (usable_w - inner_gap * (n - 1)) / n
        box_h = LAYER_H - 0.50
        for j, (title, caption, marker) in enumerate(boxes):
            x = LAYER_X0 + 0.25 + j * (box_w + inner_gap)
            y = y0 + 0.25
            ax.add_patch(FancyBboxPatch(
                (x, y), box_w, box_h,
                boxstyle="round,pad=0.02,rounding_size=0.08",
                facecolor=fill, edgecolor="black", linewidth=0.6,
                alpha=0.95, zorder=3,
            ))
            ax.text(x + box_w / 2, y + box_h - 0.22, title,
                    ha="center", va="top",
                    fontsize=11.5, fontweight="bold", color=txt, zorder=4)
            ax.text(x + box_w / 2, y + box_h - 0.65, caption,
                    ha="center", va="top",
                    fontsize=9, color=txt, zorder=4)
            if marker == "mock":
                ax.add_patch(FancyBboxPatch(
                    (x + box_w - 0.62, y + box_h - 0.32), 0.52, 0.26,
                    boxstyle="round,pad=0.01,rounding_size=0.04",
                    facecolor="#fff7ed", edgecolor="#c2410c", linewidth=0.6, zorder=5,
                ))
                ax.text(x + box_w - 0.36, y + box_h - 0.19, "mock",
                        ha="center", va="center",
                        fontsize=7.5, fontweight="bold", color="#c2410c", zorder=6)

    # --- Down arrows + their labels (in inter-layer whitespace) ------------
    arrow_x = (LAYER_X0 + LAYER_X1) / 2
    for i in range(3):
        top = layer_y[i][0]      # bottom edge of layer i
        bot = layer_y[i + 1][1]  # top edge of layer i+1
        # Arrow with small inset from band edges
        a = FancyArrowPatch(
            (arrow_x, top - 0.08), (arrow_x, bot + 0.08),
            arrowstyle="-|>", mutation_scale=18,
            color="#374151", linewidth=1.8, zorder=7,
        )
        ax.add_patch(a)
        # Label sits to the right of the arrow, vertically centered in the gap
        ax.text(arrow_x + 0.22, (top + bot) / 2, ARROW_LABELS[i],
                fontsize=9.5, color="#374151", va="center", ha="left",
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                          edgecolor="#d1d5db", linewidth=0.5))

    # --- Feedback loop (right margin, dashed green) ------------------------
    fb_x = LAYER_X1 + 0.55
    # Pass the arrow up through all 3 inter-layer gaps so it visually loops
    # from Layer 4 (bottom) back to Layer 3 (top of the world model).
    fb_top = layer_y[3][1] + 0.05    # leaves Layer 4 (top edge)
    fb_bot = layer_y[2][0] - 0.05    # enters Layer 3 (bottom edge)
    ax.add_patch(FancyArrowPatch(
        (fb_x, fb_top), (fb_x, fb_bot),
        arrowstyle="-|>", mutation_scale=16,
        color="#15803d", linewidth=2.0,
        linestyle=(0, (5, 3)), zorder=7,
    ))
    # Label centered in the right margin, vertical, color-matched
    ax.text(W - RIGHT_MARGIN / 2 + 0.20, (fb_top + fb_bot) / 2,
            "observations → refit GP + REINFORCE on ΔHV",
            fontsize=10, color="#15803d", ha="center", va="center",
            fontweight="bold", rotation=90)

    # --- Honesty footer ----------------------------------------------------
    fig.text(0.99, 0.015,
             "v1 release · pretrained policy + cached ESM3 embeddings are pre-computed indices · "
             "wet-lab simulator is the load-bearing mock",
             ha="right", fontsize=7.5, color="#666", style="italic")

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
