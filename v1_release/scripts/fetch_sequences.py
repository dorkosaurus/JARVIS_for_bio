"""Fetch AAV2 VP1 wildtype protein sequence from NCBI."""

import sys
from pathlib import Path

from Bio import Entrez, SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config

Entrez.email = config.NCBI_EMAIL


def fetch_aav2_vp1() -> tuple[str, str]:
    handle = Entrez.efetch(
        db="protein",
        id=config.AAV2_ACCESSION,
        rettype="fasta",
        retmode="text",
    )
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq), record.description


def write_fasta(seq: str, header: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")


def main() -> None:
    seq, desc = fetch_aav2_vp1()
    header = f"AAV2_WT|accession={config.AAV2_ACCESSION}|length={len(seq)}"
    write_fasta(seq, header, config.AAV2_FASTA)
    print(f"Wrote {config.AAV2_FASTA}")
    print(f"  accession:   {config.AAV2_ACCESSION}")
    print(f"  description: {desc}")
    print(f"  length:      {len(seq)} aa")
    print(f"  first 40 aa: {seq[:40]}")
    print(f"  last 20 aa:  ...{seq[-20:]}")


if __name__ == "__main__":
    main()
