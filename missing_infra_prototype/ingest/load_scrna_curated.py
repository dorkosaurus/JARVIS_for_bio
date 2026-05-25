"""Load curated gene_programs, gene_program_members, and differential_expression.

For v0 these tables are hand-curated from published AMD scRNA-seq atlases
(Voigt 2019, Orozco 2020, Menon 2019) rather than computed end-to-end from
count matrices. Each row carries a citation in the source CSV. Swapping in
real NMF / DESeq2 / MAST pipeline outputs later is a single-script change.
"""

from __future__ import annotations

import csv
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

from .init_db import DB_PATH

PROTOTYPE_ROOT = Path(__file__).resolve().parent.parent
CURATED_DIR = PROTOTYPE_ROOT / "data" / "curated"
PROGRAMS_CSV = CURATED_DIR / "gene_programs.csv"
MEMBERS_CSV = CURATED_DIR / "gene_program_members.csv"
DE_CSV = CURATED_DIR / "differential_expression.csv"


def _read_csv(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f))


def load(db_path: Path = DB_PATH) -> dict[str, int]:
    now = datetime.now(timezone.utc).isoformat()

    programs = _read_csv(PROGRAMS_CSV)
    members = _read_csv(MEMBERS_CSV)
    de = _read_csv(DE_CSV)

    program_rows = [
        (
            r["program_id"],
            r["program_label"],
            r["cell_type_label"],
            r["cell_type_ontology"],
            r["atlas_dataset"],
            r["method"],
            float(r["coherence"]) if r.get("coherence") else None,
            f"{r['source_dataset']} :: {r['citation']}",
            now,
        )
        for r in programs
    ]

    member_rows = [
        (
            r["program_id"],
            r["gene_symbol"],
            float(r["loading"]),
            int(r["rank"]),
        )
        for r in members
    ]

    de_rows = [
        (
            r["gene_symbol"],
            r["cell_type_label"],
            r["cell_type_ontology"],
            r["comparison"],
            float(r["log2fc"]),
            float(r["padj"]),
            r["method"],
            r["atlas_dataset"],
            int(r["n_case"]) if r.get("n_case") else None,
            int(r["n_control"]) if r.get("n_control") else None,
            f"curated_v0 :: {r['citation']}",
            now,
        )
        for r in de
    ]

    with sqlite3.connect(db_path) as conn:
        # Order matters because gene_program_members FK -> gene_programs.
        conn.execute("DELETE FROM gene_program_members")
        conn.execute("DELETE FROM gene_programs")
        conn.execute("DELETE FROM differential_expression")
        conn.executemany(
            "INSERT INTO gene_programs "
            "(program_id, program_label, cell_type_label, cell_type_ontology, "
            " atlas_dataset, method, coherence, source_dataset, retrieved_at) "
            "VALUES (?,?,?,?,?,?,?,?,?)",
            program_rows,
        )
        conn.executemany(
            "INSERT INTO gene_program_members "
            "(program_id, gene_symbol, loading, rank) VALUES (?,?,?,?)",
            member_rows,
        )
        conn.executemany(
            "INSERT INTO differential_expression "
            "(gene_symbol, cell_type_label, cell_type_ontology, comparison, "
            " log2fc, padj, method, atlas_dataset, n_case, n_control, "
            " source_dataset, retrieved_at) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            de_rows,
        )
        conn.commit()

    return {
        "programs": len(program_rows),
        "program_members": len(member_rows),
        "de_rows": len(de_rows),
    }


if __name__ == "__main__":
    stats = load()
    for k, v in stats.items():
        print(f"  {k}: {v}")
