# AMD hypothesis-generation: end-to-end demo run

This is the canonical demo run of the JARVIS prototype: a fresh Claude Code
session is given an AMD GWAS objective and walks the 7-step DAG by calling
MCP tools only. **No raw count matrices, GWAS sumstats, coloc inputs, or
ontology dumps are ever in the agent's context.** Every fact in the final
hypothesis is backed by a tool-returned provenance string.

The transcript below is reproduced from a live run against
`prototype/indices/jarvis.sqlite` (built from the ingest scripts under
`prototype/ingest/`). Tool inputs and outputs are abbreviated for legibility
but are unmodified in shape.

---

## Setup

`.mcp.json` at the repo root exposes two stdio servers:

- `jarvis-indices` — 8 query tools over the 5 indices (credible sets, coloc,
  pathway, gene programs, DE).
- `jarvis-workflow` — DAG retrieval (`list_workflows`, `get_workflow`,
  `find_workflow_for_objective`).

Start Claude Code from `/home/ubuntu/JARVIS_for_bio` and confirm both MCP
servers are listed in `/mcp`.

## Objective handed to the agent

> Using the GWAS indices, formulate a cellular-level hypothesis that links
> AMD risk variants to disease biology in the retinal pigment epithelium.
> Show your reasoning chain with provenance.

---

## Step −1 · Agent looks up the workflow

The agent's first move is *not* a database query — it asks what workflow
exists for this objective so it doesn't have to invent the analytic flow.

**Tool call** `workflow.find_workflow_for_objective(objective_substring="macular degeneration")`

```json
{
  "count": 1,
  "workflows": [{
    "workflow_id": "Run_GWAS_Discovery",
    "objective": "Analyze a GWAS for age-related macular degeneration (AMD) and formulate cellular hypotheses linking risk variants to disease biology.",
    "description": "Canonical GWAS->hypothesis flow as described in Hypothesis_Generation.md. ..."
  }]
}
```

**Tool call** `workflow.get_workflow(workflow_id="Run_GWAS_Discovery")`

Returns the 7-step DAG. The agent now knows which tools to call in which
order:

```
[0] Discover available GWAS                                tool=list_gwas_studies
[1] Retrieve credible variants                             tool=query_credible_sets
[2] Assign causal genes via QTL colocalization             tool=query_coloc
[3] Map causal genes to pathways                           tool=query_pathway_membership
[4] Find co-expression programs in the relevant cell type  tool=query_gene_programs
[5] Retrieve case/control differential expression          tool=query_differential_expression
[6] Synthesize hypothesis with full evidence chain         (agent reasoning, no tool call)
```

---

## Step 0 · Discover AMD GWAS studies

**Tool call** `list_gwas_studies(trait_substring="macular")`

| gwas_id              | n_loci | n_credible_variants | trait                                                |
|----------------------|--------|---------------------|------------------------------------------------------|
| `FINNGEN_R12_H7_AMD` | 36     | 527                 | Age-related macular degeneration (wet or dry)        |
| `GCST003219`         | 29     | 312                 | Advanced AMD (Fritsche et al. 2016)                  |
| `GCST90480048`       | 40     | 135                 | Senile macular degeneration NOS (PheCode 362.29)     |

**Agent decision:** anchor on `GCST003219` — Fritsche et al. 2016 is the
field-canonical reference dataset and has the deepest fine-mapping coverage.

---

## Step 1 · Credible variants in the anchor GWAS

**Tool call** `query_credible_sets(gwas_id="GCST003219", pip_threshold=0.5, limit=15)`

Returns 14 high-confidence credible variants (provenance:
`Open Targets Platform v26.03 credibleSets`). Highlights:

| variant_id     | locus    | PIP   | chrom |
|----------------|----------|-------|-------|
| `rs10922109`   | **CFH**  | 0.984 | 1     |
| `rs3750846`    | **ARMS2**| 0.992 | 10    |
| `rs2230199`    | **C3**   | 0.964 | 19    |
| `rs10033900`   | CFI      | 0.979 | 4     |
| `rs429358`     | APOE     | 0.970 | 19    |
| `rs140647181`  | COL8A1   | 0.998 | 3     |
| `rs13278062`   | TNFRSF10A| 0.947 | 8     |
| `rs2740488`    | ABCA1    | 0.546 | 9     |

**Agent decision:** anchor downstream on three loci with classical AMD
biology — CFH (`rs10922109`), the ARMS2/HTRA1 region (`rs3750846`), and C3
(`rs2230199`). Note `rs10922109` matches Fritsche 2016's published CFH lead
variant exactly — sanity check that the index reflects the source paper.

---

## Step 2 · Causal-gene assignment via QTL colocalization

**Tool call** `query_coloc(variants=["rs3750846","rs10922109","rs2230199"], pp4_threshold=0.9)`

50 coloc rows; provenance: `Open Targets Platform v26.03 colocalisation
(eQTL Catalogue + GTEx + UKB-PPP)`.

**At `rs10922109` (CFH locus):** 14 coloc rows. Every high-PP4 row is to
*CFH itself* — across left ventricle, skin, esophagus, stomach. No retina or
RPE eQTL is present (a real data gap, not a bug); the assignment to CFH is
nonetheless unambiguous because it is the only colocalized gene.

**At `rs3750846` (ARMS2 locus):** 33 coloc rows split between *ARMS2*
(testis eQTL, PP4≈0.998), *HTRA1-AS1* (DLPFC, fibroblast eQTL, PP4≈0.999),
and *HTRA1*. This **reproduces the field's well-known causal-gene ambiguity
at this locus** — the indices surface the ambiguity rather than hiding it.
Agent will carry both ARMS2 and HTRA1 forward as candidate causal genes.

**At `rs2230199` (C3 locus):** 3 coloc rows — CFP and VSIG4 (plasma pQTLs,
PP4≈0.998) plus a monocyte eQTL. CFP (Complement Factor Properdin)
colocalizing on the C3 lead variant is itself biologically informative.

**Causal gene set forwarded:** `{ARMS2, CFH, CFP, HTRA1, HTRA1-AS1, VSIG4}`.

---

## Step 3 · Map causal genes to pathways

**Tool call** `query_pathway_membership(gene_symbol="CFH")`

| pathway_id       | n_members | pathway_name                       |
|------------------|-----------|------------------------------------|
| `R-HSA-977606`   | 128       | **Regulation of Complement cascade** |
| `R-HSA-166658`   | 146       | Complement cascade                 |
| `R-HSA-168249`   | 1198      | Innate Immune System               |
| `R-HSA-168256`   | 2295      | Immune System                      |

Provenance: `Reactome v96 GMT release`.

**Tool call** `query_pathway_membership(gene_symbol="C3")` → 27 pathways
including the very specific *Alternative complement activation*
(`R-HSA-173736`, n=5) and *Activation of C3 and C5* (`R-HSA-174577`, n=7).

**Agent decision:** Both anchor genes (CFH and C3) converge on the
*Regulation of Complement cascade* pathway — the smallest non-trivial
pathway covering both. The umbrella `Innate Immune System` term is too
generic to be informative.

**Tool call** `query_pathway_members(pathway_id="R-HSA-977606")` → 128
members, includes the full C1q complex (C1QA/B/C), MAC components (C5–C9),
factor B/I/H, all CFHRs, regulators (CD46, CD55, CD59), and clusterin
(CLU — itself an AMD risk locus).

---

## Step 4 · Find a co-expression program in RPE (CL:0002586)

**Tool call** `query_gene_programs(genes=<55 candidate genes>, cell_type_ontology="CL:0002586", min_overlap=2)`

| program_id                    | overlap | coh  | overlapping_genes                                     |
|-------------------------------|---------|------|-------------------------------------------------------|
| `RETINA_RPE_NMF_COMPLEMENT`   | **13/13** | 0.86 | C1QA, C1QB, C1QC, C2, C3, CD46, CD55, CD59, CFB, CFH, CFHR1, CFHR3, CFI |

Provenance: `curated v0 from published retinal scRNA atlases (Voigt 2019,
Orozco 2020); see source_dataset for citations`.

**Interpretation:** *Every member* of the RPE complement co-expression
program is also a member of the pathway that contains the causal genes.
This is the bridge between germline genetic risk (steps 1–3) and cellular
state (step 5).

**Tool call** `query_gene_program_members(program_id="RETINA_RPE_NMF_COMPLEMENT")`
returns the program members ranked by NMF loading: CFH(0.42), C3(0.38),
CFB(0.31), CFI(0.29), C1QA(0.24)… — i.e. CFH and C3, the two anchor causal
genes, are the top two drivers of the program.

---

## Step 5 · Case/control differential expression in RPE

**Tool call** `query_differential_expression(genes=<program∪causal>, cell_type_ontology="CL:0002586", comparison="AMD_vs_control")`

11 rows, padj ≤ 0.05; summary `{mean_log2fc: +1.15, n_up: 10, n_down: 1}`.

| gene    | dir | log2FC | padj    |
|---------|-----|--------|---------|
| C3      | ↑   | +1.85  | 9.0e-05 |
| HTRA1   | ↑   | +1.78  | 1.8e-04 |
| C1QA    | ↑   | +1.65  | 2.1e-04 |
| C1QB    | ↑   | +1.50  | 3.5e-04 |
| C1QC    | ↑   | +1.48  | 3.8e-04 |
| CFH     | ↑   | +1.42  | 3.0e-04 |
| CFB     | ↑   | +1.21  | 1.7e-03 |
| CFI     | ↑   | +0.93  | 4.2e-03 |
| CFHR1   | ↑   | +0.71  | 1.2e-02 |
| CFHR3   | ↑   | +0.62  | 2.4e-02 |
| CD46    | ↓   | -0.48  | 3.8e-02 |

Provenance: `curated_v0 :: Orozco LD et al. 2020 Cell Rep 30:1246`.

**Interpretation:** The complement program is coordinately upregulated in
AMD RPE, with the strongest effects in the activators (C3, C1q, CFB,
HTRA1) and a *down*regulation in the membrane regulator CD46 — a
mechanistically coherent signature of complement dysregulation in the
disease-relevant cell type.

---

## Step 6 · Synthesized hypothesis (agent reasoning, no tool call)

> **Hypothesis.** Common AMD risk variants at CFH (`rs10922109`,
> PIP=0.984) and C3 (`rs2230199`, PIP=0.964) act through a single,
> coherent disease mechanism: dysregulated activation of the alternative
> complement cascade in the retinal pigment epithelium. The CFH risk
> allele reduces inhibitory regulation; the C3 risk allele alters the
> activatable substrate; together they shift the local complement
> equilibrium toward activation. The ARMS2/HTRA1 locus (`rs3750846`,
> PIP=0.992) contributes a parallel insult — HTRA1 upregulation in RPE
> (log2FC=+1.78, padj=1.8e-04) is concordant with the eQTL evidence and
> with HTRA1's known role in Bruch's-membrane proteolysis.
>
> **Cellular-level prediction.** RPE cells in AMD donors will show
> coordinated upregulation of a 13-gene complement co-expression program
> (top drivers CFH, C3, CFB, CFI, C1Q complex) together with
> down-regulation of the inhibitory membrane regulator CD46 — a
> mechanistically coherent activation signature.
>
> **Evidence chain.**
> 1. Variant → causal gene: `rs10922109` → CFH (coloc PP4=0.994,
>    OT v26.03); `rs2230199` → CFP/C3 region (PP4=0.998, plasma pQTL);
>    `rs3750846` → ARMS2/HTRA1-AS1 (PP4≥0.997). Causal-gene ambiguity at
>    the ARMS2 locus is acknowledged.
> 2. Gene → pathway: CFH ∈ *Regulation of Complement cascade*
>    (R-HSA-977606, n=128, Reactome v96). C3 ∈ *Alternative complement
>    activation* (R-HSA-173736, n=5, Reactome v96).
> 3. Pathway → cellular program: 13/13 members of `RETINA_RPE_NMF_COMPLEMENT`
>    are pathway members (curated v0 :: Voigt 2019, Orozco 2020).
> 4. Cellular program → disease state: 10/11 members significantly
>    upregulated in AMD RPE (mean log2FC=+1.15; curated v0 :: Orozco
>    2020 Cell Rep 30:1246).
>
> **Caveats surfaced by the indices.**
> - No retina/RPE/choroid eQTL is available in OT v26.03 for any of the
>   anchor loci; causal-gene assignment for CFH depends on a tissue-wide
>   consistent eQTL effect rather than a tissue-matched signal.
> - The ARMS2 vs HTRA1 question cannot be resolved from these indices
>   alone; both are plausible from the colocalisation and HTRA1 has
>   stronger DE support in RPE, but a definitive assignment requires a
>   RPE-matched QTL not yet in the index.

---

## What this run demonstrates

1. **Context discipline.** Every fact the hypothesis cites came back as a
   small structured row from an index. No raw GWAS sumstats, no coloc
   input matrices, no count matrices, no ontology dumps entered the
   conversation.
2. **Provenance throughout.** Each tool result carries a `source_dataset`
   or `provenance` string; the agent can — and does — cite these
   verbatim in the final hypothesis.
3. **Real biology emerges from real data.** The CFH lead variant matches
   the Fritsche 2016 published lead; the ARMS2/HTRA1 ambiguity is
   reproduced rather than papered over; the complement program in RPE
   is coordinately upregulated in AMD donors.
4. **The workflow tells the agent what to do.** The agent retrieves the
   DAG once at the start and follows it; it does not have to invent the
   analytic flow or remember which index to query when.

## How to re-run

```bash
source ~/venv/bin/activate

# (one-time) rebuild the indices from source
python -m prototype.ingest.init_db
python -m prototype.ingest.load_pathway_membership
python -m prototype.ingest.load_credible_sets_and_coloc
python -m prototype.ingest.load_scrna_curated
python -m prototype.ingest.load_workflows

# verify the indices answer correctly without MCP transport
python -m prototype.mcp_servers._smoke_test_indices

# launch Claude Code from the repo root; both MCP servers come up via .mcp.json
cd /home/ubuntu/JARVIS_for_bio
claude
```
