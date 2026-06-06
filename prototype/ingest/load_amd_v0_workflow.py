"""Seed the workflows table with Run_AMD_Hypothesis_v0.

The v0 substack-launch workflow: hands the agent a GWAS sumstats file and walks
it through Open Targets L2G + ESM3 variant viz + PaperClip literature +
v0 mock DE/pathway, ending in an agent-composed per-target mechanistic hypothesis.

Run via:
    ~/venv/bin/python -m prototype.ingest.load_amd_v0_workflow
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path

from .init_db import DB_PATH

WORKFLOW_ID = "Run_AMD_Hypothesis_v0"

WORKFLOW = {
    "workflow_id": WORKFLOW_ID,
    "objective": (
        "Generate mechanistic per-target hypotheses from an AMD GWAS using real "
        "Open Targets + ESM3 + PaperClip — fully MCP-driven"
    ),
    "description": (
        "v0 end-to-end demo: GWAS sumstats -> OT credible sets -> OT L2G gene "
        "assignment -> ESM3 fold + variant viz -> PaperClip literature -> per-target "
        "mechanistic hypothesis reports. DE & pathway steps remain honest mocks "
        "(jarvis-indices) for v0; real backends in v1."
    ),
}

STEPS = [
    {
        "step_id": "step_0_read_gwas",
        "step_order": 0,
        "step_name": "Read the GWAS sumstats file",
        "description": (
            "Open the *.sumstats.tsv file the user provided; extract study_id from "
            "## headers and the top-loci preview. Do NOT re-run fine-mapping; the "
            "OT credible sets are pre-computed."
        ),
        "tool_name": "(local file read)",
        "input_schema": {"sumstats_path": "string"},
        "output_schema": {
            "study_id": "GCST*",
            "top_loci_preview": [
                {"snp": "rsID", "chr": "int", "bp": "int", "locus_label": "str"}
            ],
        },
        "depends_on": "",
    },
    {
        "step_id": "step_1_study",
        "step_order": 1,
        "step_name": "Look up the study in Open Targets",
        "description": (
            "Confirm the GWAS exists in OT's study catalog. Returns trait, ancestry, "
            "sample size, PMID — anchors the analysis."
        ),
        "tool_name": "jarvis-ot.study_lookup",
        "input_schema": {"study_id": "GCST*"},
        "output_schema": {
            "found": "bool",
            "study": {
                "traitFromSource": "str",
                "nCases": "int",
                "nControls": "int",
                "pubmedId": "str",
            },
        },
        "depends_on": "step_0_read_gwas",
    },
    {
        "step_id": "step_2_credible_sets",
        "step_order": 2,
        "step_name": "Retrieve fine-mapped credible sets",
        "description": (
            "Pull all SuSiE credible sets for the study. Each studyLocusId becomes "
            "the key that joins variants to genes in step 3."
        ),
        "tool_name": "jarvis-ot.credible_sets_for_study",
        "input_schema": {"study_id": "GCST*"},
        "output_schema": {
            "credible_sets": [
                {"studyLocusId": "hash", "chromosome": "str", "position": "int"}
            ]
        },
        "depends_on": "step_1_study",
    },
    {
        "step_id": "step_3_l2g",
        "step_order": 3,
        "step_name": "Locus-to-gene assignment",
        "description": (
            "Get top L2G predictions joined with gene symbols. L2G score blends 29 "
            "features (eQTL/pQTL/sQTL coloc, distance, VEP, gene density). Top 6 "
            "genes by score become the target set."
        ),
        "tool_name": "jarvis-ot.l2g_top_genes",
        "input_schema": {
            "study_id": "GCST*",
            "score_threshold": "float",
            "limit": "int",
        },
        "output_schema": {
            "top_genes": [
                {
                    "studyLocusId": "hash",
                    "gene_symbol": "str",
                    "l2g_score": "float",
                }
            ]
        },
        "depends_on": "step_2_credible_sets",
    },
    {
        "step_id": "step_4_metadata_and_variants",
        "step_order": 4,
        "step_name": "Per-target: gene metadata + lead variant",
        "description": (
            "For each top L2G gene, look up canonical UniProt + transcript via "
            "jarvis-ot.gene_metadata, then jarvis-ot.lead_variant_for_locus to get "
            "the highest-PIP variant (the proximal causal candidate)."
        ),
        "tool_name": "jarvis-ot.gene_metadata + jarvis-ot.lead_variant_for_locus",
        "input_schema": {"symbol": "str", "studyLocusId": "hash"},
        "output_schema": {
            "gene": {
                "ensg": "str",
                "uniprot": "str",
                "canonical_transcript": "str",
            },
            "lead_variant": {
                "variantId": "str",
                "posteriorProbability": "float",
            },
        },
        "depends_on": "step_3_l2g",
    },
    {
        "step_id": "step_5_esm3",
        "step_order": 5,
        "step_name": "Per-target: ESM3 fold + functional consequence + PNG render",
        "description": (
            "jarvis-esm3.score_target orchestrates: VEP variant consequence -> ESM3 "
            "Forge fold (PDB + pLDDT + pTM + InterPro annotations) -> PyMOL headless "
            "render (variant residue highlighted if coding, full protein otherwise). "
            "Cached per UniProt; ~25-30s per fresh target."
        ),
        "tool_name": "jarvis-esm3.score_target",
        "input_schema": {
            "gene_symbol": "str",
            "uniprot_id": "str",
            "variant_id": "str?",
        },
        "output_schema": {
            "variant_consequence": "object",
            "fold": {"mean_plddt": "float", "pTM": "float"},
            "render": {"png_path": "str"},
        },
        "depends_on": "step_4_metadata_and_variants",
    },
    {
        "step_id": "step_6_de_mock",
        "step_order": 6,
        "step_name": "Per-target: case/control differential expression (v0 mock)",
        "description": (
            "v0 honest stub: jarvis-indices.query_differential_expression returns "
            "mock case/control DE in retina/RPE/choroid cell types. Real backends "
            "in v1."
        ),
        "tool_name": "jarvis-indices.query_differential_expression",
        "input_schema": {"gene_symbol": "str"},
        "output_schema": {
            "rows": [
                {"cell_type": "str", "log2fc": "float", "padj": "float"}
            ]
        },
        "depends_on": "step_5_esm3",
    },
    {
        "step_id": "step_7_pathway_mock",
        "step_order": 7,
        "step_name": "Per-target: pathway membership (v0 mock)",
        "description": (
            "v0 honest stub: jarvis-indices.query_pathway_membership returns "
            "Reactome/Wiki pathway memberships. Real backends in v1."
        ),
        "tool_name": "jarvis-indices.query_pathway_membership",
        "input_schema": {"gene_symbol": "str"},
        "output_schema": {
            "rows": [{"pathway_id": "str", "pathway_name": "str"}]
        },
        "depends_on": "step_5_esm3",
    },
    {
        "step_id": "step_8_literature",
        "step_order": 8,
        "step_name": "Per-target: literature corroboration",
        "description": (
            "jarvis-paperclip.literature_for_gene with disease_context='age-related "
            "macular degeneration'. Returns top N papers (PMC/bioRxiv/medRxiv/arXiv) "
            "with title, authors, paper_id, summary."
        ),
        "tool_name": "jarvis-paperclip.literature_for_gene",
        "input_schema": {
            "gene_symbol": "str",
            "disease_context": "str",
            "n": "int",
        },
        "output_schema": {
            "papers": [
                {
                    "title": "str",
                    "paper_id": "str",
                    "url": "str",
                    "summary": "str",
                }
            ]
        },
        "depends_on": "step_5_esm3",
    },
    {
        "step_id": "step_9_compose",
        "step_order": 9,
        "step_name": "Per-target: compose mechanistic hypothesis report",
        "description": (
            "Reasoning step (agent, no tool). Integrate L2G score + top features, "
            "variant consequence + PNG, mock DE + pathway notes, top papers. Write "
            "to output/<symbol>_<ensg>.md with cited provenance from every prior "
            "step. This is the ONLY step that requires reasoning; everything else "
            "is retrieval."
        ),
        "tool_name": "(agent reasoning, no tool call)",
        "input_schema": {"all_prior_step_outputs": "object"},
        "output_schema": {
            "report_path": "str",
            "evidence_chain": ["str"],
            "provenance": ["str"],
        },
        "depends_on": "step_6_de_mock,step_7_pathway_mock,step_8_literature",
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
