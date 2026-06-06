"""JARVIS Open Targets MCP server.

Backed by DuckDB over the Open Targets Platform parquets at
/home/ubuntu/JARVIS_for_bio/data/ot/. Single-connection, read-only.

Run via:
    ~/venv/bin/python -m prototype.mcp_servers.ot_server
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import duckdb
from mcp.server.fastmcp import FastMCP

DATA = Path("/home/ubuntu/JARVIS_for_bio/data/ot")
CACHE_DB = DATA / "cache.duckdb"
OT_RELEASE = "2026-03"
PROVENANCE = f"Open Targets Platform release {OT_RELEASE}"

mcp = FastMCP("jarvis-ot")

_con: Optional[duckdb.DuckDBPyConnection] = None


def _db() -> duckdb.DuckDBPyConnection:
    """Lazy-attached DuckDB connection.

    Hot-path queries (study lookup, top L2G genes, gene symbol lookup) hit
    the pre-joined slim tables in CACHE_DB. Struct-heavy detail (credible-set
    `locus` member arrays, L2G `features` SHAP arrays, `target.transcripts`)
    is left on parquet and reached via on-demand views — these are only
    called once per gene-of-interest, not per query.
    """
    global _con
    if _con is None:
        _con = duckdb.connect(":memory:")
        # Slim hot-path tables: indexed, materialized, fast
        _con.execute(f"ATTACH '{CACHE_DB}' AS ot (READ_ONLY)")
        _con.execute("CREATE VIEW study AS SELECT * FROM ot.study")
        _con.execute("CREATE VIEW target_slim AS SELECT * FROM ot.target")
        _con.execute("CREATE VIEW study_l2g_slim AS SELECT * FROM ot.study_l2g_slim")
        # Struct-heavy detail tables: parquet-backed, scanned on demand
        _con.execute(
            f"CREATE VIEW credible_set_locus AS "
            f"SELECT studyLocusId, chromosome, position, variantId, region, locus "
            f"FROM read_parquet('{DATA}/credible_set/*.parquet')"
        )
        _con.execute(
            f"CREATE VIEW l2g_features AS "
            f"SELECT studyLocusId, geneId, score, shapBaseValue, features "
            f"FROM read_parquet('{DATA}/l2g_prediction/*.parquet')"
        )
        _con.execute(
            f"CREATE VIEW target_full AS "
            f"SELECT id, approvedSymbol, approvedName, biotype, "
            f"       canonicalTranscript, genomicLocation, transcripts "
            f"FROM read_parquet('{DATA}/target/*.parquet')"
        )
    return _con


def _rows(cur) -> list[dict[str, Any]]:
    cols = [d[0] for d in cur.description]
    return [dict(zip(cols, r)) for r in cur.fetchall()]


@mcp.tool()
def study_lookup(study_id: str) -> dict[str, Any]:
    """Look up a GWAS study by its accession (e.g. GCST003219).

    Returns trait, PMID, sample size, ancestry, and disease mappings.
    """
    cur = _db().execute(
        """
        SELECT studyId, traitFromSource,
               pubmedId, publicationFirstAuthor, publicationDate,
               publicationTitle, nCases, nControls, nSamples,
               studyType, hasSumstats
        FROM study WHERE studyId = ?
        """,
        [study_id],
    )
    rows = _rows(cur)
    if not rows:
        return {"found": False, "study_id": study_id, "provenance": PROVENANCE}
    return {"found": True, "study": rows[0], "provenance": PROVENANCE}


@mcp.tool()
def credible_sets_for_study(study_id: str) -> dict[str, Any]:
    """List all fine-mapped credible sets for a GWAS study.

    Returns one row per studyLocusId with chromosome, position. Hot-path
    served from the slim cache (no struct columns).
    """
    cur = _db().execute(
        """
        SELECT DISTINCT studyLocusId, chromosome, position
        FROM study_l2g_slim WHERE studyId = ?
        ORDER BY chromosome, position
        """,
        [study_id],
    )
    rows = _rows(cur)
    return {
        "study_id": study_id,
        "count": len(rows),
        "credible_sets": rows,
        "provenance": f"{PROVENANCE} credible_set (SuSiE fine-mapping)",
    }


@mcp.tool()
def l2g_top_genes(
    study_id: str,
    score_threshold: float = 0.5,
    limit: int = 20,
) -> dict[str, Any]:
    """Top L2G (locus-to-gene) predictions for a GWAS study, joined with gene symbols.

    Each row: studyLocusId, ENSG, gene symbol, L2G score, and locus position.
    L2G score blends coloc (eQTL/pQTL/sQTL), distance, VEP, and gene density features.
    """
    cur = _db().execute(
        """
        SELECT studyLocusId, chromosome, position, geneId, gene_symbol, biotype,
               l2g_score
        FROM study_l2g_slim
        WHERE studyId = ? AND l2g_score >= ?
        ORDER BY l2g_score DESC
        LIMIT ?
        """,
        [study_id, score_threshold, limit],
    )
    rows = _rows(cur)
    return {
        "study_id": study_id,
        "score_threshold": score_threshold,
        "count": len(rows),
        "top_genes": rows,
        "provenance": f"{PROVENANCE} l2g_prediction (29-feature gradient-boosted classifier)",
    }


@mcp.tool()
def l2g_feature_contributions(
    study_locus_id: str,
    gene_id: str,
) -> dict[str, Any]:
    """Return the 29 L2G feature contributions for a (studyLocusId, geneId) pair.

    Shows what evidence the L2G model used (eQTL coloc, distance, VEP, etc).
    SHAP-style contributions sum to (score - shapBaseValue).
    """
    cur = _db().execute(
        """
        SELECT studyLocusId, geneId, score, shapBaseValue, features
        FROM l2g_features
        WHERE studyLocusId = ? AND geneId = ?
        """,
        [study_locus_id, gene_id],
    )
    rows = _rows(cur)
    if not rows:
        return {"found": False, "provenance": PROVENANCE}
    r = rows[0]
    feats = [
        {"name": f["name"], "value": f["value"], "shap": f["shapValue"]}
        for f in (r.get("features") or [])
    ]
    feats.sort(key=lambda x: abs(x["shap"] or 0), reverse=True)
    return {
        "found": True,
        "studyLocusId": r["studyLocusId"],
        "geneId": r["geneId"],
        "score": r["score"],
        "shapBaseValue": r["shapBaseValue"],
        "top_contributing_features": feats[:10],
        "all_features": feats,
        "provenance": f"{PROVENANCE} l2g_prediction features (SHAP contributions)",
    }


@mcp.tool()
def lead_variant_for_locus(study_locus_id: str) -> dict[str, Any]:
    """Return the highest-PIP variant in a credible set, with genomic coords.

    Pulls from the locus[] struct inside credible_set. Includes is95CredibleSet
    membership and posterior probability.
    """
    cur = _db().execute(
        """
        SELECT studyLocusId, chromosome, position, variantId, region, locus
        FROM credible_set_locus WHERE studyLocusId = ?
        """,
        [study_locus_id],
    )
    rows = _rows(cur)
    if not rows:
        return {"found": False, "studyLocusId": study_locus_id, "provenance": PROVENANCE}
    r = rows[0]
    members = r.get("locus") or []
    summarised = []
    for m in members:
        summarised.append({
            "variantId": m.get("variantId"),
            "posteriorProbability": m.get("posteriorProbability"),
            "is95CredibleSet": m.get("is95CredibleSet"),
            "is99CredibleSet": m.get("is99CredibleSet"),
            "logBF": m.get("logBF"),
            "beta": m.get("beta"),
            "standardError": m.get("standardError"),
        })
    summarised.sort(
        key=lambda x: x["posteriorProbability"] if x["posteriorProbability"] is not None else -1,
        reverse=True,
    )
    lead = summarised[0] if summarised else None
    return {
        "found": True,
        "studyLocusId": r["studyLocusId"],
        "chromosome": r["chromosome"],
        "position": r["position"],
        "lead_variant": lead,
        "credible_set_size": len(summarised),
        "credible_set_top5": summarised[:5],
        "provenance": f"{PROVENANCE} credible_set.locus struct (SuSiE PIPs)",
    }


@mcp.tool()
def gene_metadata(symbol_or_ensg: str) -> dict[str, Any]:
    """Look up a gene's metadata: ENSG, symbol, biotype, canonical UniProt,
    canonical transcript, chromosome / position, and approved name.

    Accepts either an HGNC symbol (e.g. "C3") or ENSG ID (e.g. "ENSG00000125730").
    """
    is_ensg = symbol_or_ensg.startswith("ENSG")
    where = "id = ?" if is_ensg else "approvedSymbol = ?"
    cur = _db().execute(
        f"""
        SELECT id, approvedSymbol, approvedName, biotype,
               canonicalTranscript, genomicLocation, transcripts
        FROM target_full WHERE {where}
        """,
        [symbol_or_ensg],
    )
    rows = _rows(cur)
    if not rows:
        return {"found": False, "query": symbol_or_ensg, "provenance": PROVENANCE}
    r = rows[0]
    canon = r.get("canonicalTranscript") or {}
    transcripts = r.get("transcripts") or []
    uniprot = None
    for t in transcripts:
        if t.get("transcriptId") == canon.get("id"):
            uniprot = t.get("uniprotId")
            break
    if uniprot is None and transcripts:
        for t in transcripts:
            if t.get("uniprotId"):
                uniprot = t.get("uniprotId")
                break
    loc = r.get("genomicLocation") or {}
    return {
        "found": True,
        "ensg": r["id"],
        "symbol": r["approvedSymbol"],
        "name": r["approvedName"],
        "biotype": r["biotype"],
        "uniprot": uniprot,
        "canonical_transcript_id": canon.get("id"),
        "chromosome": loc.get("chromosome"),
        "start": loc.get("start"),
        "end": loc.get("end"),
        "strand": loc.get("strand"),
        "provenance": f"{PROVENANCE} target (Ensembl + UniProt)",
    }


if __name__ == "__main__":
    mcp.run()
