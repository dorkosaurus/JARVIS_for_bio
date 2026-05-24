"""Smoke-test for indices_server tools, run synchronously without the MCP transport.

Calls each registered tool with realistic args and prints a one-line summary.
This catches SQL bugs and shape problems before the MCP stdio handshake.
"""

from __future__ import annotations

import asyncio
import json
from prototype.mcp_servers import indices_server as ix


async def main() -> None:
    # FastMCP wraps the original function in a Tool; .call_tool returns
    # (content_list, structured_dict). The underlying functions are also
    # still callable directly via the module — easier to test.

    print("== list_gwas_studies(trait_substring='macular') ==")
    r = ix.list_gwas_studies(trait_substring="macular")
    print(f"  {r['count']} studies")
    for s in r["studies"][:3]:
        print(f"    {s['gwas_id']:<25} loci={s['n_loci']:<3} variants={s['n_credible_variants']}")

    print("\n== query_credible_sets(locus='CFH', pip_threshold=0.5) ==")
    r = ix.query_credible_sets(locus="CFH", pip_threshold=0.5)
    print(f"  {r['count']} rows; top: {r['rows'][0]['variant_id']} PIP={r['rows'][0]['pip']:.3f}")

    print("\n== query_coloc(variants=['rs10922109'], pp4_threshold=0.9) ==")
    r = ix.query_coloc(variants=["rs10922109"], pp4_threshold=0.9)
    print(f"  {r['count']} coloc rows")
    for row in r["rows"][:5]:
        print(f"    {row['gene_symbol']:<8} {row['qtl_type']:<5} {row['tissue_label']:<35} PP4={row['pp4']:.3f}")

    print("\n== query_pathway_membership('CFH') ==")
    r = ix.query_pathway_membership("CFH")
    print(f"  {r['count']} pathways")
    for row in r["rows"]:
        print(f"    {row['pathway_id']:<18} n={row['n_members']:<4} {row['pathway_name']}")

    print("\n== query_pathway_members('R-HSA-166658') first 10 ==")
    r = ix.query_pathway_members("R-HSA-166658", limit=10)
    print(f"  total={r['count']} sample={r['members'][:10]}")

    print("\n== query_gene_programs(genes=['CFH','C3','CFB','CFI'], cell_type_ontology='CL:0002586') ==")
    r = ix.query_gene_programs(genes=["CFH", "C3", "CFB", "CFI"], cell_type_ontology="CL:0002586")
    print(f"  {r['count']} programs")
    for row in r["rows"]:
        print(f"    {row['program_id']:<35} overlap={row['n_overlap']}/{row['member_count']} coh={row['coherence']} -> {row['overlapping_genes']}")

    print("\n== query_differential_expression(genes=['CFH','C3','CFB','CFI','C1QA','C1QB','C1QC'], cell_type_ontology='CL:0002586') ==")
    r = ix.query_differential_expression(
        genes=["CFH", "C3", "CFB", "CFI", "C1QA", "C1QB", "C1QC"],
        cell_type_ontology="CL:0002586",
    )
    print(f"  {r['count']} DE rows; summary={r['summary']}")
    for row in r["rows"]:
        print(f"    {row['gene_symbol']:<8} log2FC={row['log2fc']:+.2f}  padj={row['padj']:.2e}")


if __name__ == "__main__":
    asyncio.run(main())
