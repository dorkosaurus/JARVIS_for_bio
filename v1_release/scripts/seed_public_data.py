"""Seed the public pre-training CSVs used to fit the GP world model.

The four CSVs are synthetic-but-grounded: target values are hand-authored
against published anchors (AAV2 ~0.015 transduction / 0.10 escape /
0.15 inflammation; AAV.7m8 ~0.45 / 0.30 / 0.35; etc.) and the VP1
sequences are generated to match the mutation classes the papers describe
(insertion variants for Dalkara directed-evolution intermediates;
VR-IV/VR-VIII substitutions for Byrne 4D-R100 family; etc.).

This is the *public data layer*. v2+ would replace each CSV with the
actual paper's published readouts where machine-readable.
"""

import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config


# ----------------------------- helpers -----------------------------------


def load_anchor_sequences() -> tuple[str, str]:
    aav2 = str(SeqIO.read(config.AAV2_FASTA, "fasta").seq)
    aav7m8 = str(SeqIO.read(config.AAV7M8_FASTA, "fasta").seq)
    return aav2, aav7m8


def apply_spec(spec: tuple, aav2_seq: str, aav7m8_seq: str) -> tuple[str, dict]:
    """Return (vp1_sequence, engineered_feature_dict)."""
    kind = spec[0]
    if kind == "anchor":
        if spec[1] == "AAV2":
            return aav2_seq, {
                "has_7mer_insertion": False,
                "insertion_length": 0,
                "insertion_peptide": None,
                "hamming_to_aav2": 0,
            }
        if spec[1] == "AAV7m8":
            return aav7m8_seq, {
                "has_7mer_insertion": True,
                "insertion_length": len(config.AAV7M8_PEPTIDE),
                "insertion_peptide": config.AAV7M8_PEPTIDE,
                "hamming_to_aav2": len(config.AAV7M8_PEPTIDE),
            }
        raise ValueError(f"unknown anchor: {spec[1]}")
    if kind == "insertion":
        peptide = spec[1]
        pos = config.INSERTION_AFTER_RESIDUE
        new_seq = aav2_seq[:pos] + peptide + aav2_seq[pos:]
        return new_seq, {
            "has_7mer_insertion": True,
            "insertion_length": len(peptide),
            "insertion_peptide": peptide,
            "hamming_to_aav2": len(peptide),
        }
    if kind == "substitution":
        muts = spec[1]  # list[(pos_1indexed, new_aa)]
        new = list(aav2_seq)
        for pos, aa in muts:
            new[pos - 1] = aa
        return "".join(new), {
            "has_7mer_insertion": False,
            "insertion_length": 0,
            "insertion_peptide": None,
            "hamming_to_aav2": len(muts),
        }
    raise ValueError(f"unknown spec kind: {kind}")


# ----------------------------- data ---------------------------------------

# Each row: (capsid_id, mutation_spec, rpe_transduction, neut_escape,
#           inflammation_score, extras_dict, source_label)
# Anchors (AAV2 / AAV.7m8) appear in multiple CSVs with the same capsid_id
# so the embedding index can dedupe.

DALKARA_2013 = [
    ("AAV2",   ("anchor", "AAV2"),    0.015, 0.10, 0.15),
    ("AAV7m8", ("anchor", "AAV7m8"),  0.45,  0.30, 0.35),
    ("DAL_001", ("insertion", "GRPSDSP"),   0.05, 0.12, 0.17),
    ("DAL_002", ("insertion", "LRPSDSP"),   0.07, 0.15, 0.19),
    ("DAL_003", ("insertion", "AGRPSDP"),   0.10, 0.17, 0.22),
    ("DAL_004", ("insertion", "GETRPSP"),   0.15, 0.20, 0.25),
    ("DAL_005", ("insertion", "LRGETRP"),   0.20, 0.22, 0.27),
    ("DAL_006", ("insertion", "LAGETRP"),   0.28, 0.25, 0.30),
    ("DAL_007", ("insertion", "ALGETRP"),   0.32, 0.26, 0.31),
    ("DAL_008", ("insertion", "DGETTRP"),   0.25, 0.23, 0.28),
    ("DAL_009", ("insertion", "LALGGTRP"),  0.40, 0.29, 0.34),
    ("DAL_010", ("insertion", "LAGETRPS"),  0.35, 0.27, 0.32),
    ("DAL_011", ("insertion", "ALGETTPR"),  0.30, 0.26, 0.30),
    ("DAL_012", ("insertion", "LALGGRP"),   0.36, 0.27, 0.32),
    ("DAL_013", ("insertion", "LAGGETRP"),  0.38, 0.28, 0.33),
]

BYRNE_2020 = [
    ("AAV2",   ("anchor", "AAV2"),    0.015, 0.10, 0.15),
    ("AAV7m8", ("anchor", "AAV7m8"),  0.45,  0.30, 0.35),
    # 4D-R100 family: VR-IV / VR-VIII substitutions
    ("BYR_001", ("substitution", [(449, "V"), (454, "R"), (588, "K")]), 0.22, 0.32, 0.24),
    ("BYR_002", ("substitution", [(461, "K"), (588, "K"), (590, "T")]), 0.28, 0.38, 0.28),
    ("BYR_003", ("substitution", [(449, "S"), (457, "N"), (461, "K")]), 0.18, 0.28, 0.22),
    ("BYR_004", ("substitution", [(454, "R"), (461, "K"), (467, "S"), (590, "T")]), 0.35, 0.42, 0.30),
    ("BYR_005", ("substitution", [(449, "V"), (461, "K"), (588, "K"), (590, "T")]), 0.40, 0.48, 0.34),
    ("BYR_006", ("substitution", [(449, "V"), (461, "K"), (584, "N")]), 0.25, 0.35, 0.26),
    ("BYR_007", ("substitution", [(457, "N"), (588, "K")]), 0.20, 0.30, 0.22),
    ("BYR_008", ("substitution", [(449, "V"), (454, "R"), (461, "K"), (467, "S"), (590, "T")]), 0.48, 0.52, 0.38),
    ("BYR_009", ("substitution", [(454, "R"), (461, "K"), (590, "T")]), 0.30, 0.40, 0.28),
    ("BYR_010", ("substitution", [(449, "V"), (461, "K"), (467, "S"), (588, "K")]), 0.38, 0.45, 0.32),
]

KOTTERMAN_2015 = [
    # Only neut_escape and n_serum_samples are reported.
    ("AAV2",    ("anchor", "AAV2"),                      0.10, 200),
    ("AAV7m8",  ("anchor", "AAV7m8"),                    0.30, 200),
    ("KOT_001", ("insertion", "RGDSTPR"),                0.40, 180),
    ("KOT_002", ("insertion", "NNPAHQD"),                0.55, 200),
    ("KOT_003", ("substitution", [(449, "T"), (467, "S"), (588, "K"), (708, "K")]), 0.45, 200),
    ("KOT_004", ("substitution", [(454, "R"), (461, "K"), (467, "S")]), 0.25, 180),
    ("KOT_005", ("insertion", "QSRTAGR"),                0.65, 150),
    ("KOT_006", ("substitution", [(265, "T"), (449, "V")]), 0.15, 200),
]

REICHEL_BUCHER = [
    # Only inflammation_score and ifn_beta_fold are reported.
    ("AAV2",   ("anchor", "AAV2"),                       0.15, 1.0),
    ("AAV7m8", ("anchor", "AAV7m8"),                     0.35, 4.2),
    ("REI_001", ("insertion", "LALGETTRP"),              0.40, 5.1),  # ADVM-022 derivative
    ("REI_002", ("substitution", [(449, "V"), (461, "K"), (467, "S")]), 0.28, 3.0),  # 4D-R100 family
    ("REI_003", ("substitution", [(454, "R"), (588, "K")]), 0.22, 2.4),
    ("REI_004", ("insertion", "RAATRPN"),                0.45, 6.0),
    ("REI_005", ("substitution", [(449, "V"), (454, "R"), (461, "K"), (467, "S"), (505, "M"), (590, "T")]), 0.55, 9.5),
    ("REI_006", ("insertion", "GRPSDSPA"),               0.30, 3.6),
    ("REI_007", ("substitution", [(265, "T"), (705, "K")]), 0.10, 1.2),
    ("REI_008", ("substitution", [(449, "V"), (454, "R"), (461, "K"), (467, "S"), (584, "N"), (588, "K"), (708, "K")]), 0.60, 12.0),
]


# ----------------------------- builders -----------------------------------


def build_dalkara(aav2: str, aav7m8: str) -> tuple[pd.DataFrame, dict]:
    rows = []
    seqs = {}
    for cid, spec, t, e, i in DALKARA_2013:
        seq, eng = apply_spec(spec, aav2, aav7m8)
        seqs[cid] = seq
        rows.append(
            {
                "capsid_id": cid,
                "vp1_sequence": seq,
                "rpe_transduction": t,
                "neut_escape": e,
                "inflammation_score": i,
                **eng,
                "source": "Dalkara 2013",
                "source_doi": "10.1126/scitranslmed.3005708",
            }
        )
    return pd.DataFrame(rows), seqs


def build_byrne(aav2: str, aav7m8: str) -> tuple[pd.DataFrame, dict]:
    rows = []
    seqs = {}
    for cid, spec, t, e, i in BYRNE_2020:
        seq, eng = apply_spec(spec, aav2, aav7m8)
        seqs[cid] = seq
        rows.append(
            {
                "capsid_id": cid,
                "vp1_sequence": seq,
                "rpe_transduction": t,
                "neut_escape": e,
                "inflammation_score": i,
                **eng,
                "source": "Byrne 2020",
                "source_doi": "10.1172/JCI133380",
            }
        )
    return pd.DataFrame(rows), seqs


def build_kotterman(aav2: str, aav7m8: str) -> tuple[pd.DataFrame, dict]:
    rows = []
    seqs = {}
    for cid, spec, e, n in KOTTERMAN_2015:
        seq, eng = apply_spec(spec, aav2, aav7m8)
        seqs[cid] = seq
        rows.append(
            {
                "capsid_id": cid,
                "vp1_sequence": seq,
                "neut_escape": e,
                "n_serum_samples": n,
                **eng,
                "source": "Kotterman 2015",
                "source_doi": "10.1038/gt.2014.115",
            }
        )
    return pd.DataFrame(rows), seqs


def build_reichel_bucher(aav2: str, aav7m8: str) -> tuple[pd.DataFrame, dict]:
    rows = []
    seqs = {}
    for cid, spec, i, ifn in REICHEL_BUCHER:
        seq, eng = apply_spec(spec, aav2, aav7m8)
        seqs[cid] = seq
        rows.append(
            {
                "capsid_id": cid,
                "vp1_sequence": seq,
                "inflammation_score": i,
                "ifn_beta_fold": ifn,
                **eng,
                "source": "Reichel 2017 / Bucher 2021",
                "source_doi": "10.1016/j.ymthe.2017.07.010",
            }
        )
    return pd.DataFrame(rows), seqs


def write_public_fasta(unique_seqs: dict[str, str], out_path: Path) -> None:
    """One FASTA record per unique capsid_id.

    Sequence-level dedup is handled downstream by the embedder's hash cache.
    Multiple capsid_ids with the same sequence (e.g. REI_001's LALGETTRP
    insertion equals AAV.7m8) all need their own index row pointing at the
    shared embedding."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        for cid, seq in unique_seqs.items():
            f.write(f">{cid}|public_data\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def main() -> None:
    aav2, aav7m8 = load_anchor_sequences()
    config.PUBLIC_DATA_DIR.mkdir(parents=True, exist_ok=True)

    builders = [
        ("dalkara_2013.csv", build_dalkara),
        ("byrne_2020.csv", build_byrne),
        ("kotterman_2015.csv", build_kotterman),
        ("reichel_bucher_inflammation.csv", build_reichel_bucher),
    ]

    all_seqs: dict[str, str] = {}
    total_rows = 0
    for fname, builder in builders:
        df, seqs = builder(aav2, aav7m8)
        df.to_csv(config.PUBLIC_DATA_DIR / fname, index=False)
        all_seqs.update(seqs)
        total_rows += len(df)
        print(f"  {fname}: {len(df)} rows")

    public_fasta = config.PUBLIC_DATA_DIR / "public_sequences.fasta"
    write_public_fasta(all_seqs, public_fasta)

    n_unique = len({s for s in all_seqs.values()})
    print(f"\n  total rows across CSVs: {total_rows}")
    print(f"  unique capsid_ids:      {len(all_seqs)}")
    print(f"  unique sequences:       {n_unique}")
    print(f"  public FASTA:           {public_fasta}")


if __name__ == "__main__":
    main()
