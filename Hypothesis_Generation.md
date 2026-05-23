# AI-driven hypothesis generation

## Purpose

This document describes how we can concretely enable agentic driven hypothesis generation in the context of biological discovery.

---

## The user experience we want

```
User:  I have a GWAS for AMD, can you do your best to analyze this dataset and formulate 
       hypotheses for how these findings affect AMD disease biology at a cellular level?

Agent: Looking
       [...thinking...]  # 2 minutes pass

Agent: From this GWAS, I was able to find that variant rs1234 that colocalizes with a 
       trans-pQTL for gene ABC which is a part of the WNT signalling pathway. I was also 
       able to find single cell evidence that suggests this gene and many others in the WNT 
       signalling pathway are co-expressed and upregulated in RPE cells (retinal pigment 
       epithelial cells) in patients with AMD vs those that do not. Given this, I was able 
       to form the hypothesis that upregulation of gene ABC leads AMD via upregulation of 
       the WNT signalling pathway in RPE cells.
```

---

## How it's done today

GWAS → Hypothesis:

1. Create credible sets of variants within the GWAS.
2. Perform causal gene analysis (including pQTL identification).
3. Identify pathways it plays a role in — these become gene sets we want to interrogate.
4. Identify gene programs from single-cell expression data where we can find these gene programs.
5. Perform differential expression between case/controls to determine if pathways are differentially expressed between people who have AMD and people who do not.
6. Form the hypothesis.

### Why an agent can't perform these tasks today

**Step 1: Credible set generation.** An agent handed a raw GWAS summary statistic file needs to run fine-mapping software (SuSiE, FINEMAP) to generate credible sets. These tools require compute infrastructure, reference LD panels matched to the study population, and parameter choices that affect downstream results. No agent can invoke these reliably without a pre-configured execution environment and validated reference data. The agent doesn't fail because it doesn't know what to do — it fails because the tools aren't callable.

**Step 2: Causal gene assignment and pQTL identification.** Colocalization against a pQTL dataset requires knowing which pQTL dataset to use, which tissues are relevant, and what statistical thresholds are appropriate. The eQTL Catalogue contains hundreds of datasets. An agent querying it naively will either retrieve too much (token economics barrier) or make arbitrary choices about which dataset is relevant without the context to justify them.

**Step 3: Pathway identification.** Mapping a gene to pathways is tractable — databases like WikiPathways and Reactome are queryable. But the output is a gene set, and the agent now needs to decide which gene sets are biologically meaningful in the context of AMD specifically, not just which ones contain gene ABC. That judgment requires context the agent doesn't have unless it's been pre-loaded.

**Step 4: Gene program identification from single-cell data.** This is where the token economics barrier becomes insurmountable without pre-computed infrastructure. Identifying gene programs from a single-cell atlas requires access to count matrices, dimensionality reduction, and gene program analysis (NMF, topic modeling). The relevant atlas for AMD — cells from retinal tissue across AMD and control donors — may be hundreds of gigabytes. An agent cannot retrieve and reason over this at query time.

**Step 5: Differential expression.** Running differential expression between AMD and control RPE cells requires the same count matrices, a statistical framework (DESeq2, edgeR), and choices about covariates, normalization, and multiple testing correction. Not callable without pre-configured infrastructure and pre-selected datasets.

**Step 6: Hypothesis formation.** This is the one step an agent can do well — if it has been given the outputs of steps 1–5 in a compact, well-structured form. The hypothesis formation step is a reasoning task. Every step before it is a data engineering task. Agents are being asked to do both, and the data engineering is what breaks them.

---

### An agent-native architecture to enable efficient hypothesis formation

```
Objective --> Analytic Flow -->
                            --> Analytic Step 1 --> Look up precomputed index
                                                    Or query a foundational model
                            --> Analytic Step 2 --> Look up precomputed index
                                                    Or query a foundational model
                            --> Analytic Step N --> Look up precomputed index
                                                    Or query a foundational model
```

The key architectural insight is that every step in a biological analytic workflow is a data retrieval problem, not a compute problem. The compute happened — or should have happened — at ingest time. The agent's job is to retrieve pre-computed results, reason over them, and decide what to retrieve next. This is tractable. Recomputing everything from raw data at query time is not.

**The analytic flow graph.** The agent begins with an objective: analyze a GWAS for AMD and formulate cellular hypotheses. It looks up the canonical analytic flow for this objective in a workflow graph — a DAG where each node is a well-defined task with known inputs, outputs, and data dependencies. The human genetics workflow is a known pattern. The agent doesn't need to invent the workflow. It follows it.

**Step 1: Credible set retrieval.** Rather than running fine-mapping software, the agent queries a pre-computed credible set index. The GWAS summary statistics have already been processed against a reference LD panel. The agent receives a ranked list of credible variants with posterior inclusion probabilities. Input: GWAS identifier. Output: credible variant set. Token cost: hundreds of tokens, not billions.

```
Agent action: query_credible_sets(gwas_id="AMD_2024", threshold=0.95)
Returns: [rs1234 (PIP=0.87), rs5678 (PIP=0.72), ...]
```

**Step 2: Causal gene assignment and pQTL colocalization.** The agent queries a pre-computed colocalization index — eQTL Catalogue results materialized across relevant tissues — filtered to retinal and ocular tissues by ontology-grounded retrieval. "Retinal pigment epithelial cells" resolves via Cell Ontology to every relevant tissue in the index before a single record is retrieved. The agent receives colocalization probabilities for each credible variant against known pQTLs. Input: credible variant set, tissue ontology term. Output: colocalized gene assignments with evidence provenance.

```
Agent action: query_coloc(variants=["rs1234"], tissue="CL:0002586", pp4_threshold=0.8)
Returns: rs1234 colocalizes with pQTL for gene ABC (PP4=0.91, dataset=GTEx_retina_v9)
```

**Step 3: Pathway context.** The agent queries WikiPathways and Reactome for gene ABC. It receives a pre-summarized evidence hierarchy: the pathways gene ABC participates in, the other members of those pathways, and a confidence score based on experimental evidence. The WNT signaling pathway is returned as the highest-confidence association. Token cost: a structured JSON record, not a literature scan.

```
Agent action: query_pathway_membership(gene="ABC")
Returns: WNT signaling pathway (evidence_score=0.94), VEGF signaling (0.61), ...
```

**Step 4: Gene program identification.** The agent queries a pre-computed gene program index built over a retinal single-cell atlas — NMF components and topic model outputs materialized at ingest, indexed by cell type and gene set membership. It retrieves the gene programs active in RPE cells that contain gene ABC and other WNT pathway members. No count matrices are retrieved. No dimensionality reduction is run at query time. Input: gene set, cell type ontology term. Output: gene program memberships with co-expressed genes.

```
Agent action: query_gene_programs(genes=["ABC", ...WNT_members], cell_type="CL:0002586")
Returns: Program 7 (WNT/RPE, coherence=0.88): ABC, DEF, GHI, ...
```

**Step 5: Differential expression retrieval.** The agent queries a pre-computed differential expression index — DESeq2 results materialized across AMD vs. control comparisons in retinal datasets, stored in Parquet and queryable via DuckDB. It retrieves the differential expression statistics for gene program 7 members in RPE cells. Input: gene set, cell type, case/control comparison. Output: log fold changes, adjusted p-values, dataset provenance.

```
Agent action: query_differential_expression(genes=["ABC", "DEF", "GHI"], 
              cell_type="CL:0002586", comparison="AMD_vs_control")
Returns: ABC (log2FC=1.8, padj=0.0001), DEF (log2FC=2.1, padj=0.0003), ...
```

**Step 6: Hypothesis formation.** The agent now has everything it needs: a credible variant, a colocalized gene, a pathway, a co-expressed gene program, and differential expression evidence in the relevant cell type. This is the one step that is purely a reasoning task. The agent synthesizes the evidence chain into a structured hypothesis with full provenance.

```
Agent action: Generates Hypothesis
Returns:
    Hypothesis: Upregulation of gene ABC drives AMD pathology via activation of the
    WNT signaling pathway in retinal pigment epithelial cells.

    Evidence chain:
    - rs1234 (PIP=0.87) colocalizes with pQTL for ABC (PP4=0.91, GTEx_retina_v9)
    - ABC is a high-confidence member of the WNT signaling pathway (score=0.94, WikiPathways)
    - ABC and WNT pathway members are co-expressed in RPE gene program 7 (coherence=0.88)
    - Program 7 genes are upregulated in AMD RPE cells vs controls
      (mean log2FC=1.9, padj<0.001, dataset=AMD_scRNA_atlas_2023)
```

**What made this possible.** Every intermediate result was retrieved from a pre-computed index in milliseconds. The agent never saw raw count matrices, GWAS summary statistics, or colocalization inputs. It saw structured, provenance-tagged answers to well-defined questions. The total token cost of the workflow is thousands of tokens — not trillions. The hypothesis formation step, the only step that required genuine reasoning, was given exactly the evidence it needed and nothing more.

---

## Prototyping an agent-native environment

```
Objective -> Agent (Claude Code)
                --> Workflow DB via MCP server
                --> Pre-computed indices (SQLite) via MCP server
                --> Foundational models via MCP server
```

This prototype instantiates the AMD hypothesis generation workflow using components that are either publicly available or constructible from public data. The goal is to demonstrate that the architecture works end to end, not to build production infrastructure.

### 1. The agentic interface: Claude Code

Claude Code acts as the orchestrating agent. It receives the user objective, looks up the analytic flow, and issues tool calls against the MCP servers at each step. It never touches raw data directly.

### 2. Pre-computed indices (SQLite via MCP)

Each step in the AMD workflow maps to a SQLite table, constructible from public sources:

| Table | Source | Content |
|---|---|---|
| credible_sets | GWAS Catalog + SuSiE pre-run | variant, gwas_id, PIP, locus |
| coloc_results | eQTL Catalogue + coloc pre-run | variant, gene, tissue_ontology, PP4, dataset |
| pathway_membership | WikiPathways + Reactome | gene, pathway, evidence_score |
| gene_programs | retinal scRNA atlas + NMF pre-run | program_id, cell_type_ontology, genes, coherence |
| differential_expression | AMD scRNA atlas + DESeq2 pre-run | gene, cell_type_ontology, comparison, log2FC, padj, dataset |

Each table is exposed as a read-only MCP tool. Claude Code generates the query, the MCP server executes it, and the result comes back as a structured JSON record. Token cost per step: hundreds of tokens.

### 3. Foundational models (via MCP)

For steps where embedding-based retrieval is preferable to structured lookup — gene program identification being the primary case — a foundational model endpoint is exposed via MCP. In v0 this is a lightweight wrapper around a pre-computed embedding index (cell type embeddings from a public retinal atlas), queryable by gene set similarity.

### 4. Workflow DB (SQLite via MCP)

A single table encodes the AMD analytic flow as a DAG — nodes are steps, edges are dependencies, each node carries its tool name, input schema, and output schema. Claude Code reads this at the start of each session to know what workflow to follow and in what order.

```sql
workflow_steps: step_id, step_name, tool_name, input_schema, output_schema, depends_on
```

### What this prototype will demonstrate

A complete end-to-end run of the AMD hypothesis generation workflow — from GWAS identifier to structured hypothesis with full evidence provenance — using only pre-computed indices and a frontier model for reasoning. If it works, it validates the architecture. If it breaks, it tells us exactly which component needs to be hardened for production.
