"""JARVIS indices MCP server.

Exposes the five pre-computed indices (credible_sets, coloc_results,
pathway_membership, gene_programs, differential_expression) as MCP tools.
Each tool maps to a parametrised SQL query and returns structured JSON
with provenance fields so the agent can build a hypothesis with citations.

Run via:
    ~/venv/bin/python -m prototype.mcp_servers.indices_server
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, Optional

from mcp.server.fastmcp import FastMCP

DB_PATH = (
    Path(__file__).resolve().parent.parent / "indices" / "jarvis.sqlite"
)

mcp = FastMCP("jarvis-indices")


def _conn() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def _rows(cur: sqlite3.Cursor) -> list[dict[str, Any]]:
    return [dict(r) for r in cur.fetchall()]


# ---------------------------------------------------------------------------
# Step 0: Discovery — let the agent find what's actually in the indices
# ---------------------------------------------------------------------------

@mcp.tool()
def list_gwas_studies(trait_substring: Optional[str] = None) -> dict[str, Any]:
    """List GWAS studies present in the credible_sets index.

    Args:
        trait_substring: optional case-insensitive substring filter on trait.
    Returns:
        {studies: [{gwas_id, trait, n_credible_variants, n_loci}], count}
    """
    sql = (
        "SELECT gwas_id, trait, COUNT(*) AS n_credible_variants, "
        "       COUNT(DISTINCT locus) AS n_loci "
        "FROM credible_sets "
    )
    args: list[Any] = []
    if trait_substring:
        sql += "WHERE LOWER(trait) LIKE ? "
        args.append(f"%{trait_substring.lower()}%")
    sql += "GROUP BY gwas_id, trait ORDER BY n_credible_variants DESC"
    with _conn() as c:
        rows = _rows(c.execute(sql, args))
    return {"count": len(rows), "studies": rows}


# ---------------------------------------------------------------------------
# Step 1: credible_sets
# ---------------------------------------------------------------------------

@mcp.tool()
def query_credible_sets(
    gwas_id: Optional[str] = None,
    locus: Optional[str] = None,
    pip_threshold: float = 0.5,
    limit: int = 50,
) -> dict[str, Any]:
    """Retrieve fine-mapped credible variants.

    Args:
        gwas_id: filter to one GWAS study (e.g. "GCST003219").
        locus: filter to a gene-anchored locus label (e.g. "CFH").
        pip_threshold: minimum posterior inclusion probability.
        limit: max rows to return.
    Returns:
        {count, rows: [{variant_id, gwas_id, locus, pip, chrom, pos, ref, alt,
                        fine_mapper, trait, source_dataset}], provenance}
    """
    sql = "SELECT * FROM credible_sets WHERE pip >= ? "
    args: list[Any] = [pip_threshold]
    if gwas_id:
        sql += "AND gwas_id = ? "
        args.append(gwas_id)
    if locus:
        sql += "AND locus = ? "
        args.append(locus)
    sql += "ORDER BY pip DESC LIMIT ?"
    args.append(limit)
    with _conn() as c:
        rows = _rows(c.execute(sql, args))
    return {
        "count": len(rows),
        "rows": rows,
        "provenance": "Open Targets Platform v26.03 credibleSets",
    }


# ---------------------------------------------------------------------------
# Step 2: coloc_results
# ---------------------------------------------------------------------------

@mcp.tool()
def query_coloc(
    variants: Optional[list[str]] = None,
    gene_symbol: Optional[str] = None,
    tissue_substring: Optional[str] = None,
    qtl_type: Optional[str] = None,
    pp4_threshold: float = 0.8,
    limit: int = 100,
) -> dict[str, Any]:
    """Retrieve molecular-QTL colocalization for credible variants.

    Args:
        variants: list of rsIDs to look up.
        gene_symbol: filter to a specific colocalized gene.
        tissue_substring: case-insensitive substring on tissue_label.
        qtl_type: one of "eqtl", "pqtl", "sqtl", "sceqtl".
        pp4_threshold: minimum coloc PP4 (P(shared causal variant)).
        limit: max rows to return.
    Returns:
        {count, rows, provenance}
    """
    sql = "SELECT * FROM coloc_results WHERE pp4 >= ? "
    args: list[Any] = [pp4_threshold]
    if variants:
        placeholders = ",".join("?" * len(variants))
        sql += f"AND variant_id IN ({placeholders}) "
        args.extend(variants)
    if gene_symbol:
        sql += "AND gene_symbol = ? "
        args.append(gene_symbol)
    if tissue_substring:
        sql += "AND LOWER(tissue_label) LIKE ? "
        args.append(f"%{tissue_substring.lower()}%")
    if qtl_type:
        sql += "AND qtl_type = ? "
        args.append(qtl_type)
    sql += "ORDER BY pp4 DESC LIMIT ?"
    args.append(limit)
    with _conn() as c:
        rows = _rows(c.execute(sql, args))
    return {
        "count": len(rows),
        "rows": rows,
        "provenance": "Open Targets Platform v26.03 colocalisation (eQTL Catalogue + GTEx + UKB-PPP)",
    }


# ---------------------------------------------------------------------------
# Step 3: pathway membership
# ---------------------------------------------------------------------------

@mcp.tool()
def query_pathway_membership(gene_symbol: str) -> dict[str, Any]:
    """Retrieve curated pathways a gene belongs to.

    Args:
        gene_symbol: HGNC symbol, e.g. "CFH".
    Returns:
        {count, rows: [{pathway_id, pathway_name, pathway_db, evidence_score, n_members, ...}], provenance}
    """
    sql = (
        "SELECT pathway_id, pathway_name, pathway_db, evidence_score, "
        "       n_members, source_dataset "
        "FROM pathway_membership WHERE gene_symbol = ? "
        "ORDER BY evidence_score DESC, n_members ASC"
    )
    with _conn() as c:
        rows = _rows(c.execute(sql, [gene_symbol]))
    return {
        "count": len(rows),
        "rows": rows,
        "provenance": "Reactome v96 GMT release",
    }


@mcp.tool()
def query_pathway_members(pathway_id: str, limit: int = 200) -> dict[str, Any]:
    """List gene members of a pathway.

    Args:
        pathway_id: Reactome stable ID, e.g. "R-HSA-166658".
    Returns:
        {pathway_id, count, members: [gene_symbol, ...]}
    """
    with _conn() as c:
        rows = _rows(
            c.execute(
                "SELECT gene_symbol FROM pathway_members "
                "WHERE pathway_id = ? ORDER BY gene_symbol LIMIT ?",
                [pathway_id, limit],
            )
        )
    return {
        "pathway_id": pathway_id,
        "count": len(rows),
        "members": [r["gene_symbol"] for r in rows],
        "provenance": "Reactome v96 GMT release",
    }


# ---------------------------------------------------------------------------
# Step 4: gene programs
# ---------------------------------------------------------------------------

@mcp.tool()
def query_gene_programs(
    genes: Optional[list[str]] = None,
    cell_type_ontology: Optional[str] = None,
    min_overlap: int = 1,
) -> dict[str, Any]:
    """Find scRNA gene programs that contain a set of genes in a given cell type.

    Args:
        genes: gene symbols to look for as program members.
        cell_type_ontology: Cell Ontology ID, e.g. "CL:0002586" (RPE).
        min_overlap: minimum number of input genes that must appear in the program.
    Returns:
        {count, rows: [{program_id, program_label, cell_type_label, cell_type_ontology,
                        atlas_dataset, coherence, n_member_overlap, overlapping_genes,
                        all_member_count}], provenance}
    """
    # SQL params bind positionally; JOIN-clause params come before WHERE,
    # which come before HAVING. Build args in that exact order.
    join_args: list[Any] = []
    where_clauses: list[str] = []
    where_args: list[Any] = []

    join_filter = ""
    if genes:
        placeholders = ",".join("?" * len(genes))
        join_filter = f"AND gpm.gene_symbol IN ({placeholders})"
        join_args.extend(genes)

    if cell_type_ontology:
        where_clauses.append("gp.cell_type_ontology = ?")
        where_args.append(cell_type_ontology)

    sql = f"""
        SELECT gp.program_id, gp.program_label, gp.cell_type_label,
               gp.cell_type_ontology, gp.atlas_dataset, gp.method, gp.coherence,
               gp.source_dataset,
               COUNT(DISTINCT CASE WHEN gpm.gene_symbol IS NOT NULL THEN gpm.gene_symbol END) AS n_overlap,
               GROUP_CONCAT(DISTINCT gpm.gene_symbol) AS overlapping_genes,
               (SELECT COUNT(*) FROM gene_program_members m WHERE m.program_id = gp.program_id) AS member_count
        FROM gene_programs gp
        LEFT JOIN gene_program_members gpm
          ON gpm.program_id = gp.program_id {join_filter}
        {('WHERE ' + ' AND '.join(where_clauses)) if where_clauses else ''}
        GROUP BY gp.program_id
        HAVING n_overlap >= ?
        ORDER BY n_overlap DESC, gp.coherence DESC
    """
    args = join_args + where_args + [min_overlap]
    with _conn() as c:
        rows = _rows(c.execute(sql, args))
    return {
        "count": len(rows),
        "rows": rows,
        "provenance": "curated v0 from published retinal scRNA atlases (Voigt 2019, Orozco 2020); see source_dataset for citations",
    }


@mcp.tool()
def query_gene_program_members(program_id: str) -> dict[str, Any]:
    """Get the ranked gene members of a single program."""
    with _conn() as c:
        rows = _rows(
            c.execute(
                "SELECT gene_symbol, loading, rank FROM gene_program_members "
                "WHERE program_id = ? ORDER BY rank",
                [program_id],
            )
        )
    return {"program_id": program_id, "count": len(rows), "members": rows}


# ---------------------------------------------------------------------------
# Step 5: differential expression
# ---------------------------------------------------------------------------

@mcp.tool()
def query_differential_expression(
    genes: list[str],
    cell_type_ontology: Optional[str] = None,
    comparison: str = "AMD_vs_control",
    padj_threshold: float = 0.05,
) -> dict[str, Any]:
    """Retrieve differential expression statistics for a gene set in a cell type.

    Args:
        genes: gene symbols to look up.
        cell_type_ontology: Cell Ontology ID, e.g. "CL:0002586" (RPE).
        comparison: case/control comparison label, default "AMD_vs_control".
        padj_threshold: only return rows with padj <= this.
    Returns:
        {count, rows, summary: {mean_log2fc, n_up, n_down}, provenance}
    """
    if not genes:
        return {"count": 0, "rows": [], "summary": {}, "provenance": ""}
    placeholders = ",".join("?" * len(genes))
    sql = (
        f"SELECT * FROM differential_expression "
        f"WHERE gene_symbol IN ({placeholders}) AND comparison = ? AND padj <= ? "
    )
    args: list[Any] = list(genes) + [comparison, padj_threshold]
    if cell_type_ontology:
        sql += "AND cell_type_ontology = ? "
        args.append(cell_type_ontology)
    sql += "ORDER BY padj"
    with _conn() as c:
        rows = _rows(c.execute(sql, args))
    log2fcs = [r["log2fc"] for r in rows]
    summary = {
        "mean_log2fc": sum(log2fcs) / len(log2fcs) if log2fcs else None,
        "n_up": sum(1 for x in log2fcs if x > 0),
        "n_down": sum(1 for x in log2fcs if x < 0),
    }
    return {
        "count": len(rows),
        "rows": rows,
        "summary": summary,
        "provenance": "curated v0 from Orozco 2020 AMD scRNA atlas",
    }


if __name__ == "__main__":
    mcp.run()
