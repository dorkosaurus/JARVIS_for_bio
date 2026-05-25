"""Ingest pathway membership from the Reactome bulk GMT release.

Source: https://reactome.org/download/current/ReactomePathways.gmt.zip
Format: one pathway per line, TAB-separated: pathway_name, R-HSA-id, gene1, gene2, ...

Strategy: download once (cached), filter to pathways that contain at least one
of our anchor AMD genes, and write to pathway_membership + pathway_members.
This keeps the index small (we only need pathways relevant to the demo loci)
while preserving full Reactome provenance.
"""

from __future__ import annotations

import sqlite3
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import requests

from .init_db import DB_PATH

PROTOTYPE_ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = PROTOTYPE_ROOT / "data" / "raw"
GMT_URL = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
GMT_LOCAL_ZIP = RAW_DIR / "ReactomePathways.gmt.zip"
GMT_MEMBER = "ReactomePathways.gmt"

# Reactome release we're pinning provenance against.
# (Manually bumped when we refresh; the API exposes the current version
# at https://reactome.org/ContentService/data/database/version)
REACTOME_RELEASE = "v96"
SOURCE_DATASET = f"Reactome_{REACTOME_RELEASE}_GMT"

# Anchor AMD-relevant genes. Any Reactome pathway containing at least one of
# these is pulled into the index. Members of those pathways are then all
# included so the agent can query for "what else is in this pathway".
ANCHOR_GENES = {
    # Complement cascade (CFH locus story)
    "CFH", "CFHR1", "CFHR2", "CFHR3", "CFHR4", "CFHR5",
    "CFI", "CFB", "CFD", "CFP",
    "C2", "C3", "C4A", "C4B", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
    "C1QA", "C1QB", "C1QC", "C1R", "C1S",
    "CD46", "CD55", "CD59",
    # ARMS2/HTRA1 locus story
    "ARMS2", "HTRA1",
    # Lipid/cholesterol (APOE locus, also AMD-associated)
    "APOE", "CETP", "LIPC", "ABCA1",
    # Extracellular matrix / Bruch's membrane
    "TIMP3", "COL8A1", "COL4A3",
    # VEGF axis (wet AMD relevance)
    "VEGFA",
}


def download_gmt(force: bool = False) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    if GMT_LOCAL_ZIP.exists() and not force:
        return GMT_LOCAL_ZIP
    resp = requests.get(GMT_URL, timeout=60)
    resp.raise_for_status()
    GMT_LOCAL_ZIP.write_bytes(resp.content)
    return GMT_LOCAL_ZIP


def parse_gmt(zip_path: Path) -> list[tuple[str, str, list[str]]]:
    """Returns [(pathway_name, pathway_id, [gene_symbols]), ...]."""
    with zipfile.ZipFile(zip_path) as z:
        raw = z.read(GMT_MEMBER).decode("utf-8")
    pathways: list[tuple[str, str, list[str]]] = []
    for line in raw.splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        name, pid, *genes = parts
        # GMT may have duplicate gene entries within a line; dedupe but preserve order.
        seen: set[str] = set()
        deduped: list[str] = []
        for g in genes:
            if g and g not in seen:
                seen.add(g)
                deduped.append(g)
        pathways.append((name, pid, deduped))
    return pathways


def load(db_path: Path = DB_PATH, anchors: set[str] = ANCHOR_GENES) -> dict[str, int]:
    zip_path = download_gmt()
    all_pathways = parse_gmt(zip_path)

    # Filter to pathways that contain at least one anchor gene.
    relevant = [
        (name, pid, genes)
        for (name, pid, genes) in all_pathways
        if any(g in anchors for g in genes)
    ]

    now = datetime.now(timezone.utc).isoformat()

    membership_rows: list[tuple] = []
    members_rows: list[tuple] = []
    for name, pid, genes in relevant:
        n = len(genes)
        for g in genes:
            # We only register pathway_membership rows for anchor genes — keeps
            # the gene→pathway lookup table focused on AMD-relevant genes. The
            # full member list (pathway_members) holds everything so the agent
            # can still ask "what else is in this pathway?".
            if g in anchors:
                membership_rows.append(
                    (
                        g,
                        pid,
                        name,
                        "Reactome",
                        1.0,  # Reactome curator-assigned membership; binary in GMT.
                        n,
                        SOURCE_DATASET,
                        now,
                    )
                )
            members_rows.append((pid, g))

    with sqlite3.connect(db_path) as conn:
        conn.execute("DELETE FROM pathway_membership WHERE pathway_db = 'Reactome'")
        conn.execute(
            "DELETE FROM pathway_members WHERE pathway_id IN ("
            "SELECT DISTINCT pathway_id FROM pathway_members WHERE pathway_id LIKE 'R-HSA-%')"
        )
        conn.executemany(
            "INSERT OR REPLACE INTO pathway_membership "
            "(gene_symbol, pathway_id, pathway_name, pathway_db, evidence_score, "
            " n_members, source_dataset, retrieved_at) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            membership_rows,
        )
        conn.executemany(
            "INSERT OR REPLACE INTO pathway_members (pathway_id, gene_symbol) "
            "VALUES (?, ?)",
            members_rows,
        )
        conn.commit()

    return {
        "pathways_total_in_gmt": len(all_pathways),
        "pathways_relevant": len(relevant),
        "membership_rows": len(membership_rows),
        "member_rows": len(members_rows),
    }


if __name__ == "__main__":
    stats = load()
    for k, v in stats.items():
        print(f"  {k}: {v}")
