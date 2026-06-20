"""Build AAV.7m8 reference sequence (AAV2 + LALGETTRP insertion at 587-588).

Per Dalkara D et al. 2013, Sci Transl Med 5:189ra76.
"""

import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config


def insert_peptide(seq: str, peptide: str, after_residue_1indexed: int) -> str:
    pos = after_residue_1indexed  # 0-indexed slice point == 1-indexed residue index
    return seq[:pos] + peptide + seq[pos:]


def write_fasta(seq: str, header: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")


def main() -> None:
    if not config.AAV2_FASTA.exists():
        raise FileNotFoundError(
            f"AAV2 wildtype FASTA missing: {config.AAV2_FASTA}. "
            "Run scripts/fetch_sequences.py first."
        )
    aav2 = SeqIO.read(config.AAV2_FASTA, "fasta")
    aav2_seq = str(aav2.seq)

    pos = config.AAV7M8_INSERTION_AFTER
    aav7m8_seq = insert_peptide(aav2_seq, config.AAV7M8_PEPTIDE, pos)

    header = (
        f"AAV7m8|base=AAV2|insert={config.AAV7M8_PEPTIDE}"
        f"@{pos}-{pos+1}|length={len(aav7m8_seq)}"
    )
    write_fasta(aav7m8_seq, header, config.AAV7M8_FASTA)

    flank_before = aav2_seq[pos - 5 : pos]
    flank_after = aav2_seq[pos : pos + 5]

    print(f"Wrote {config.AAV7M8_FASTA}")
    print(f"  AAV2 length:    {len(aav2_seq)} aa")
    print(f"  AAV.7m8 length: {len(aav7m8_seq)} aa  (+ {len(config.AAV7M8_PEPTIDE)})")
    print(f"  insertion:      {config.AAV7M8_PEPTIDE} after residue {pos}")
    print(
        f"  flanks:         "
        f"...{flank_before}-[{config.AAV7M8_PEPTIDE}]-{flank_after}..."
    )


if __name__ == "__main__":
    main()
