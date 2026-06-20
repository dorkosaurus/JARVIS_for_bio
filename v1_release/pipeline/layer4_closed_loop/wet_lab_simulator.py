"""Wet-lab simulator — the load-bearing mock for v1.

Calibrated against published intravitreal AAV literature:
  - AAV2:    transduction ~0.015, escape ~0.10, inflammation ~0.15
  - AAV.7m8: transduction ~0.45,  escape ~0.30, inflammation ~0.35
    (Dalkara 2013, Kotterman 2015, Reichel/Bucher 2017/2021,
     ADVM-022 clinical dose-limiting toxicity).

Inputs are sequence-level features only. ESM3 embeddings are not used —
biology is captured at the sequence level. In v2+ this function is
replaced by the real wet lab; the closed loop is unchanged.
"""

import sys
from difflib import SequenceMatcher
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import config


def peptide_quality(peptide: str | None) -> float:
    """0..1 similarity to AAV.7m8 reference peptide LALGETTRP.

    Uses difflib SequenceMatcher ratio — captures both identity and shared
    substrings. AAV.7m8's biology depends on the GETTRP motif; this metric
    rewards any peptide that preserves it."""
    if not peptide:
        return 0.0
    return SequenceMatcher(None, peptide, config.AAV7M8_PEPTIDE).ratio()


def fitness_penalty(hamming_to_aav2: int) -> float:
    """0..1; rises past hamming=14 (capsid folding / packaging cliff).

    Cliff raised from 10 → 14 so stacked variants (insertion + 2-4 substitutions
    on non-VIII VR loops) fit comfortably. Published combination capsids
    (4D-R100 family, peptide-display variants with multiple loop edits) suggest
    AAV2 backbones can tolerate ~13-15 total structural changes before
    packaging / virion stability fails.
    """
    excess = max(0, hamming_to_aav2 - 14)
    return float(np.tanh(excess / 6.0))


def simulate_wet_lab(
    has_7mer_insertion: bool,
    insertion_peptide: str | None,
    insertion_length: int,
    hamming_to_aav2: int,
    noise_std: float = config.SIMULATOR_NOISE_STD,
    rng: np.random.Generator | None = None,
) -> dict:
    """Simulate intravitreal AAV outcomes for one capsid variant.

    Returns dict with rpe_transduction, neut_escape, inflammation_score
    (each clipped to [0, 1]) plus simulator-version metadata."""
    if rng is None:
        rng = np.random.default_rng()

    pq = peptide_quality(insertion_peptide) if has_7mer_insertion else 0.0
    fp = fitness_penalty(hamming_to_aav2)
    fitness_factor = 1.0 - 0.6 * fp

    # Number of point substitutions on the AAV2 backbone (excluding the inserted
    # peptide, whose residues are at novel positions). For insertion-only
    # variants this is 0; for stacked variants it counts the VR-loop swaps.
    sub_count = max(0, hamming_to_aav2 - insertion_length) if has_7mer_insertion else hamming_to_aav2

    # --- Axis 1: RPE transduction (intravitreal) ---
    # AAV.7m8-like insertions (pq ~ 1.0) get the full 0.40 bonus.
    transduction = (0.015 + 0.40 * pq) * fitness_factor

    # --- Axis 2: NAb escape ---
    # Three contributions, in increasing biological specificity:
    #   - peptide novelty + insertion size (insertion-only motif)
    #   - VR-loop substitutions erode existing antibody epitopes
    #   - stacked synergy: when both present, antibody recognition collapses
    #     across multiple surfaces simultaneously
    escape_from_ins = (0.10 * pq + 0.10 * (insertion_length / 9.0)) if has_7mer_insertion else 0.0
    escape_from_subs = 0.45 * (1.0 - np.exp(-sub_count / 3.0))
    if has_7mer_insertion and sub_count > 0:
        synergy = 0.15 * np.tanh(pq * sub_count / 3.0)
    else:
        synergy = 0.0
    escape_change = escape_from_ins + escape_from_subs + synergy
    neut_escape = (0.10 + escape_change) * fitness_factor

    # --- Constraint: inflammation_score ---
    # AAV2 baseline 0.15. Rises with structural change; sharp rise past cliff.
    infl_from_change = 0.10 * (1.0 - np.exp(-hamming_to_aav2 / 5.0))
    infl_from_insertion = 0.10 * pq if has_7mer_insertion else 0.0
    infl_from_cliff = 0.40 * fp
    inflammation = (
        0.15 + infl_from_change + infl_from_insertion + infl_from_cliff
    )

    # --- Noise ---
    transduction += rng.normal(0, noise_std)
    neut_escape += rng.normal(0, noise_std)
    inflammation += rng.normal(0, noise_std)

    return {
        "rpe_transduction": float(np.clip(transduction, 0, 1)),
        "neut_escape": float(np.clip(neut_escape, 0, 1)),
        "inflammation_score": float(np.clip(inflammation, 0, 1)),
        "_simulator_version": config.SIMULATOR_VERSION,
        "_is_simulated": True,
    }


# ----------------------------- sanity check --------------------------------


def _sanity_check() -> None:
    """Run all 80 generated variants + 2 anchors through the simulator
    (noise off) and report distribution statistics."""
    from Bio import SeqIO

    rng = np.random.default_rng(0)

    # Anchors
    aav2 = simulate_wet_lab(False, None, 0, 0, noise_std=0.0, rng=rng)
    aav7m8 = simulate_wet_lab(True, config.AAV7M8_PEPTIDE, 9, 9, noise_std=0.0, rng=rng)
    print("Anchors (noise off):")
    print(
        f"  AAV2:    t={aav2['rpe_transduction']:.3f}  "
        f"e={aav2['neut_escape']:.3f}  i={aav2['inflammation_score']:.3f}"
    )
    print(
        f"  AAV.7m8: t={aav7m8['rpe_transduction']:.3f}  "
        f"e={aav7m8['neut_escape']:.3f}  i={aav7m8['inflammation_score']:.3f}"
    )

    # All 80 generated variants
    rows = []
    for rec in SeqIO.parse(config.CAPSID_VARIANTS_FASTA, "fasta"):
        parts = dict(p.split("=") for p in rec.description.split("|")[1:] if "=" in p)
        has_ins = bool(int(parts["has_7mer"]))
        peptide = parts["insertion"] or None
        hamming = int(parts["hamming"])
        ins_len = len(peptide) if peptide else 0
        out = simulate_wet_lab(has_ins, peptide, ins_len, hamming, noise_std=0.0, rng=rng)
        rows.append(
            {
                "capsid_id": rec.id.split("|")[0],
                "class": parts["class"],
                "has_7mer": has_ins,
                "peptide": peptide,
                "hamming": hamming,
                **{k: v for k, v in out.items() if not k.startswith("_")},
            }
        )

    import pandas as pd

    df = pd.DataFrame(rows)
    print(f"\nAll 80 variants (noise off):")
    for cls in ["substitution", "insertion"]:
        sub = df[df["class"] == cls]
        print(f"  {cls} (n={len(sub)}):")
        for col in ["rpe_transduction", "neut_escape", "inflammation_score"]:
            print(
                f"    {col:<22} "
                f"min={sub[col].min():.3f}  max={sub[col].max():.3f}  "
                f"mean={sub[col].mean():.3f}"
            )

    constraint = config.INFLAMMATION_THRESHOLD
    n_passing = (df.inflammation_score < constraint).sum()
    print(f"\nConstraint (inflammation < {constraint}): {n_passing}/{len(df)} pass")
    top_t = df.nlargest(5, "rpe_transduction")[
        ["capsid_id", "class", "peptide", "rpe_transduction", "neut_escape", "inflammation_score"]
    ]
    print("\nTop 5 by transduction:")
    print(top_t.to_string(index=False))


if __name__ == "__main__":
    _sanity_check()
