"""Seed the workflows + workflow_steps tables with the AMD hypothesis-generation DAG.

This encodes the canonical analytic flow from the design doc:
  GWAS -> credible sets -> coloc -> pathway -> gene programs -> DE -> hypothesis

Each step records the MCP tool the agent should call, its input/output shape, and
its upstream dependencies. The agent reads this DAG at session start instead of
having to invent the workflow.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

from .init_db import DB_PATH

WORKFLOW_ID = "Run_GWAS_Discovery"

WORKFLOW = {
    "workflow_id": WORKFLOW_ID,
    "objective": (
        "Analyze a GWAS for age-related macular degeneration (AMD) and formulate "
        "cellular hypotheses linking risk variants to disease biology."
    ),
    "description": (
        "Canonical GWAS->hypothesis flow as described in Hypothesis_Generation.md. "
        "Each step retrieves from a pre-computed index; no raw count matrices, "
        "GWAS sumstats, or coloc inputs are ever in the agent's context. "
        "The final step is pure reasoning over compact, provenance-tagged evidence."
    ),
}

STEPS = [
    {
        "step_id": "step_0_discover",
        "step_order": 0,
        "step_name": "Discover available GWAS",
        "description": (
            "Find which AMD GWAS studies are present in the credible_sets index. "
            "Pick the largest / most authoritative study to anchor the analysis."
        ),
        "tool_name": "list_gwas_studies",
        "input_schema": {"trait_substring": "string?"},
        "output_schema": {"studies": [{"gwas_id": "str", "trait": "str", "n_loci": "int"}]},
        "depends_on": "",
    },
    {
        "step_id": "step_1_credible_sets",
        "step_order": 1,
        "step_name": "Retrieve credible variants",
        "description": (
            "For the chosen GWAS, retrieve fine-mapped credible variants (PIP >= 0.5). "
            "These are the candidate causal variants to chase downstream."
        ),
        "tool_name": "query_credible_sets",
        "input_schema": {"gwas_id": "string", "pip_threshold": "float"},
        "output_schema": {"rows": [{"variant_id": "rsID", "locus": "gene", "pip": "float"}]},
        "depends_on": "step_0_discover",
    },
    {
        "step_id": "step_2_coloc",
        "step_order": 2,
        "step_name": "Assign causal genes via QTL colocalization",
        "description": (
            "For each credible variant, retrieve molecular-QTL colocalization (PP4 >= 0.8). "
            "QTL types include eqtl, pqtl, sqtl, sceqtl. A high-PP4 coloc with a pQTL is "
            "the strongest causal-gene assignment."
        ),
        "tool_name": "query_coloc",
        "input_schema": {"variants": "list[rsID]", "pp4_threshold": "float"},
        "output_schema": {"rows": [{"variant_id": "rsID", "gene_symbol": "str", "tissue_label": "str", "pp4": "float"}]},
        "depends_on": "step_1_credible_sets",
    },
    {
        "step_id": "step_3_pathway",
        "step_order": 3,
        "step_name": "Map causal genes to pathways",
        "description": (
            "For each colocalized causal gene, retrieve curated pathway memberships. "
            "Smaller pathways (n_members < ~200) are more biologically specific than "
            "umbrella terms like 'Innate Immune System'."
        ),
        "tool_name": "query_pathway_membership",
        "input_schema": {"gene_symbol": "string"},
        "output_schema": {"rows": [{"pathway_id": "str", "pathway_name": "str", "n_members": "int"}]},
        "depends_on": "step_2_coloc",
    },
    {
        "step_id": "step_4_gene_programs",
        "step_order": 4,
        "step_name": "Find co-expression programs containing the causal-gene set",
        "description": (
            "Query scRNA gene programs WITHOUT prescribing a cell type. Pass the "
            "causal gene list (optionally expanded with the pathway's member list) "
            "and let the data surface which cell type(s) contain a coherent "
            "co-expression program for those genes. Each returned program carries "
            "its own cell_type_ontology — collect those for step 5. Do NOT pick a "
            "cell type yourself: the disease-relevant cell type(s) are an OUTPUT "
            "of this step, not an input. If multiple cell types light up, carry "
            "all of them forward."
        ),
        "tool_name": "query_gene_programs",
        "input_schema": {"genes": "list[str]", "cell_type_ontology": "OMIT — must be left unset"},
        "output_schema": {"rows": [{"program_id": "str", "program_label": "str", "cell_type_label": "str", "cell_type_ontology": "CL ID", "n_overlap": "int", "coherence": "float", "overlapping_genes": "csv"}]},
        "depends_on": "step_3_pathway",
    },
    {
        "step_id": "step_5_de",
        "step_order": 5,
        "step_name": "Retrieve case/control differential expression in the cell type(s) step 4 returned",
        "description": (
            "For each cell_type_ontology returned by step 4 (NOT one you pick), "
            "retrieve case/control differential expression for the program members. "
            "If step 4 returned multiple cell types, query each separately and "
            "compare the signatures. Coordinated up/down-regulation is the "
            "cellular-level evidence linking the genetic signal to disease state."
        ),
        "tool_name": "query_differential_expression",
        "input_schema": {"genes": "list[str]", "cell_type_ontology": "CL ID from step 4 output (per program)", "comparison": "string"},
        "output_schema": {"rows": [{"gene_symbol": "str", "log2fc": "float", "padj": "float", "cell_type_ontology": "CL ID"}]},
        "depends_on": "step_4_gene_programs",
    },
    {
        "step_id": "step_6_hypothesis",
        "step_order": 6,
        "step_name": "Synthesize hypothesis with full evidence chain",
        "description": (
            "Compose a structured hypothesis that links the credible variant -> "
            "colocalized gene -> pathway -> cell-type program -> case/control "
            "differential expression. Every claim must cite the provenance returned by "
            "the corresponding tool. This is the ONLY step that requires reasoning; "
            "every prior step was retrieval."
        ),
        "tool_name": "(agent reasoning, no tool call)",
        "input_schema": {"all_prior_step_outputs": "object"},
        "output_schema": {"hypothesis": "str", "evidence_chain": "list[str]", "provenance": "list[str]"},
        "depends_on": "step_5_de",
    },
]


def load(db_path: Path = DB_PATH) -> dict[str, int]:
    with sqlite3.connect(db_path) as conn:
        conn.execute("DELETE FROM workflow_steps WHERE workflow_id = ?", [WORKFLOW_ID])
        conn.execute("DELETE FROM workflows WHERE workflow_id = ?", [WORKFLOW_ID])
        conn.execute(
            "INSERT INTO workflows (workflow_id, objective, description) VALUES (?,?,?)",
            (WORKFLOW["workflow_id"], WORKFLOW["objective"], WORKFLOW["description"]),
        )
        for s in STEPS:
            conn.execute(
                "INSERT INTO workflow_steps "
                "(workflow_id, step_id, step_order, step_name, description, "
                " tool_name, input_schema, output_schema, depends_on) "
                "VALUES (?,?,?,?,?,?,?,?,?)",
                (
                    WORKFLOW_ID,
                    s["step_id"],
                    s["step_order"],
                    s["step_name"],
                    s["description"],
                    s["tool_name"],
                    json.dumps(s["input_schema"]),
                    json.dumps(s["output_schema"]),
                    s["depends_on"],
                ),
            )
        conn.commit()
    return {"workflows": 1, "steps": len(STEPS)}


if __name__ == "__main__":
    stats = load()
    for k, v in stats.items():
        print(f"  {k}: {v}")
