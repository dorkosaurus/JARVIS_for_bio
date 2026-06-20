"""Generate AAV2 VP1 capsid variants for the optimization candidate pool.

Two variant classes:
  - Substitution: 1-6 residue swaps at VR-IV, V, VIII, IX (Tseng &
    Agbandje-McKenna 2014 epitope/tropism hotspots).
  - Insertion: a peptide inserted between residues 587 and 588 (the
    AAV.7m8 trick, Dalkara 2013).

Peptide library is a small curated set inspired by AAV.7m8 and 4D-R100
family. v2+ would learn this with RFdiffusion + ProteinMPNN.
"""

import hashlib
import sys
from pathlib import Path

import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
import config

AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")
# Exclude cysteine from surface-loop substitutions (disulfide risk).
AA_NOCYS = [a for a in AA_ALPHABET if a != "C"]

# Curated 7-9mer peptide library. LALGETTRP (AAV.7m8, Dalkara 2013) anchors;
# other entries mimic motifs from the retinal-tropic peptide pools reported in
# Dalkara 2013 / Byrne 2020 (4D-R100 family).
INSERTION_LIBRARY = [
    "LALGETTRP",  # AAV.7m8 (Dalkara 2013)
    "LAGETRP",
    "ALGETRP",
    "DGETTRP",
    "LRGETRP",
    "VAGETRP",
    "LANGEGRP",
    "LAVGENRP",
    "LAKDETRP",
    "MELGETRP",
    "LAYGETRP",
    "LALSGTRP",
    "LALDQTR",
    "GETTRPS",
    "ALSETRP",
    "NETIGRP",
    "DSPTHPK",   # 4D-R100 family inspired
    "PDSTTRP",
    "GTSVNRP",
    "MGGTRSP",
    "GASRTPR",
    "NLAGETT",
    "RAGGTRP",
    "LIGETRH",
    "ANGETLP",
    "LAGGTTP",
    "LSTGRPS",
    "RAATPNS",
    "LALGETER",
    "LAGETRPG",
]


def capsid_id(seq: str) -> str:
    """Stable 8-char hash of sequence."""
    return hashlib.sha1(seq.encode()).hexdigest()[:8]


def substitution_residue_pool() -> list[int]:
    """1-indexed residue positions across substitution-eligible VR loops."""
    positions = []
    for loop_name in config.SUBSTITUTION_LOOPS:
        start, end = config.VR_LOOPS[loop_name]
        positions.extend(range(start, end + 1))
    return positions


def pick_alternative_aa(wt_aa: str, rng: np.random.Generator) -> str:
    pool = [a for a in AA_NOCYS if a != wt_aa]
    return str(rng.choice(pool))


def generate_substitution_variants(
    wt_seq: str, n: int, rng: np.random.Generator
) -> list[dict]:
    pool = substitution_residue_pool()
    variants: list[dict] = []
    seen: set[str] = set()
    attempts = 0
    while len(variants) < n:
        attempts += 1
        if attempts > 50 * n:
            raise RuntimeError(
                f"Could not generate {n} unique substitution variants "
                f"after {attempts} attempts."
            )
        k = int(
            rng.integers(
                config.MIN_SUBSTITUTIONS_PER_VARIANT,
                config.MAX_SUBSTITUTIONS_PER_VARIANT + 1,
            )
        )
        positions = sorted(
            int(p) for p in rng.choice(pool, size=k, replace=False)
        )
        new = list(wt_seq)
        mutations: list[str] = []
        for pos in positions:
            wt_aa = wt_seq[pos - 1]
            alt = pick_alternative_aa(wt_aa, rng)
            new[pos - 1] = alt
            mutations.append(f"{wt_aa}{pos}{alt}")
        new_seq = "".join(new)
        if new_seq in seen or new_seq == wt_seq:
            continue
        seen.add(new_seq)
        variants.append(
            {
                "sequence": new_seq,
                "class": "substitution",
                "mutations": ",".join(mutations),
                "has_7mer_insertion": False,
                "insertion_peptide": None,
                "insertion_length": 0,
                "hamming_to_aav2": k,
            }
        )
    return variants


def generate_insertion_variants(
    wt_seq: str, n: int, rng: np.random.Generator
) -> list[dict]:
    if len(INSERTION_LIBRARY) < n:
        raise ValueError(
            f"Insertion library has {len(INSERTION_LIBRARY)} peptides; "
            f"need {n}."
        )
    library = list(INSERTION_LIBRARY)
    rng.shuffle(library)
    peptides = library[:n]

    pos = config.INSERTION_AFTER_RESIDUE  # 1-indexed; slice point in 0-indexed string
    variants: list[dict] = []
    for peptide in peptides:
        new_seq = wt_seq[:pos] + peptide + wt_seq[pos:]
        variants.append(
            {
                "sequence": new_seq,
                "class": "insertion",
                "mutations": f"ins{pos}_{pos+1}:{peptide}",
                "has_7mer_insertion": True,
                "insertion_peptide": peptide,
                "insertion_length": len(peptide),
                # We treat the inserted residues as the hamming-distance proxy
                # since they're entirely new positions.
                "hamming_to_aav2": len(peptide),
            }
        )
    return variants


def stacked_residue_pool() -> list[int]:
    """1-indexed residue positions for stacked-variant substitutions (non-VIII loops)."""
    positions = []
    for loop_name in config.STACKED_SUBSTITUTION_LOOPS:
        start, end = config.VR_LOOPS[loop_name]
        positions.extend(range(start, end + 1))
    return positions


def generate_stacked_variants(
    wt_seq: str, n: int, rng: np.random.Generator
) -> list[dict]:
    """7-mer insertion at 587 + 2-4 substitutions at non-VIII VR loops.

    Stacking the AAV.7m8 trick (insertion drives transduction) with epitope
    erosion at distant VR loops (substitutions drive escape) is the move that
    can put a capsid into the target zone — both axes high simultaneously,
    without disrupting the insertion site itself.
    """
    pos = config.INSERTION_AFTER_RESIDUE
    sub_pool = stacked_residue_pool()
    library = list(INSERTION_LIBRARY)

    variants: list[dict] = []
    seen: set[str] = set()
    attempts = 0
    while len(variants) < n:
        attempts += 1
        if attempts > 50 * n:
            raise RuntimeError(
                f"Could not generate {n} unique stacked variants after {attempts} attempts."
            )
        peptide = str(rng.choice(library))
        k = int(rng.integers(config.MIN_STACKED_SUBSTITUTIONS,
                             config.MAX_STACKED_SUBSTITUTIONS + 1))
        positions = sorted(int(p) for p in rng.choice(sub_pool, size=k, replace=False))

        new = list(wt_seq)
        sub_mutations: list[str] = []
        for p in positions:
            wt_aa = wt_seq[p - 1]
            alt = pick_alternative_aa(wt_aa, rng)
            new[p - 1] = alt
            sub_mutations.append(f"{wt_aa}{p}{alt}")
        substituted = "".join(new)
        # Insert peptide between residues 587 and 588 of the substituted sequence
        stacked_seq = substituted[:pos] + peptide + substituted[pos:]

        if stacked_seq in seen:
            continue
        seen.add(stacked_seq)

        mutations = f"ins{pos}_{pos+1}:{peptide}+" + ",".join(sub_mutations)
        variants.append(
            {
                "sequence": stacked_seq,
                "class": "stacked",
                "mutations": mutations,
                "has_7mer_insertion": True,
                "insertion_peptide": peptide,
                "insertion_length": len(peptide),
                # Total structural change from AAV2: insertion length + substitutions
                "hamming_to_aav2": len(peptide) + k,
            }
        )
    return variants


def write_fasta(variants: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        for v in variants:
            cid = capsid_id(v["sequence"])
            header = (
                f">{cid}|class={v['class']}"
                f"|mutations={v['mutations']}"
                f"|has_7mer={int(v['has_7mer_insertion'])}"
                f"|insertion={v['insertion_peptide'] or ''}"
                f"|hamming={v['hamming_to_aav2']}"
            )
            f.write(header + "\n")
            seq = v["sequence"]
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def main() -> None:
    if not config.AAV2_FASTA.exists():
        raise FileNotFoundError(
            f"AAV2 wildtype FASTA missing: {config.AAV2_FASTA}. "
            "Run scripts/fetch_sequences.py first."
        )
    wt = SeqIO.read(config.AAV2_FASTA, "fasta")
    wt_seq = str(wt.seq)

    rng = np.random.default_rng(config.RNG_SEED)
    subs = generate_substitution_variants(
        wt_seq, config.N_SUBSTITUTION_VARIANTS, rng
    )
    inserts = generate_insertion_variants(
        wt_seq, config.N_INSERTION_VARIANTS, rng
    )
    stacked = generate_stacked_variants(
        wt_seq, config.N_STACKED_VARIANTS, rng
    )
    variants = subs + inserts + stacked

    write_fasta(variants, config.CAPSID_VARIANTS_FASTA)

    hamming = [v["hamming_to_aav2"] for v in variants]
    print(f"Wrote {config.CAPSID_VARIANTS_FASTA}")
    print(f"  substitution variants:  {len(subs)}")
    print(f"  insertion variants:     {len(inserts)}")
    print(f"  stacked variants:       {len(stacked)}  (insertion + non-VIII substitutions)")
    print(f"  total:                  {len(variants)}")
    print(f"  AAV2 reference length:  {len(wt_seq)} aa")
    print(f"  hamming-to-AAV2 range:  {min(hamming)}-{max(hamming)}")


if __name__ == "__main__":
    main()
