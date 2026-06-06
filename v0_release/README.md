# Drug discovery has a reasoning problem

*An inference-oriented architecture for agent-driven hypothesis generation, with a v0 that runs end-to-end on MCP.*

---

Biology spent a decade getting big. The infrastructure we built during that time was designed for human scientists to query, interpret, and act on data. What's changed is the emergence of heterogeneous AI systems capable of reasoning over that data at a scale no human team can match. The infrastructure wasn't built for them, and it shows. Without in silico intelligence that can actually plug into biological data at scale, we can't explore hypothesis space fully, and most of it goes unexplored.

The arithmetic is unforgiving. A serious AMD target-discovery effort might touch twenty GWAS loci. Each one needs fine-mapping, colocalization, differential expression, variant consequence annotation, pathway enrichment, and literature synthesis before you can say anything mechanistic about it. Most teams work three or four loci well and deprioritize the rest, not because the biology matters less, but because analyst time runs out. The hypotheses that don't get generated don't fail. They just never exist.

AI agents can change that throughput. They don't deprioritize locus six because locus one was more compelling. They can run the same reasoning protocol across all twenty candidates in parallel.

What throughput alone doesn't resolve is whether the hypotheses are any good. An agent may not know what genome-wide significance means, or when a colocalization PP4 of 0.6 is compelling versus marginal, or that ARMS2 has been contested in the AMD literature for fifteen years. It generates confidently regardless. At scale that moves the bottleneck from generation to triage. A scientist still has to sort through the output, and now there's more of it.

There's a second problem sitting underneath that one. When computation and reasoning are tangled together, fine-mapping running inline and the agent reasoning over the output, there's no clean separation between what the tool returned and what the agent inferred. You can't hand a report to a skeptical MD scientist and point to the number, the source, and the version. An agent you can't audit is an agent you can't calibrate, and calibration is what turns generated hypotheses into something a drug discovery team actually acts on.

Both problems point to the same architectural fix: separate reads from writes.

The writes, fine-mapping, colocalization, fold prediction, enrichment, happen once, in batch, ahead of hypothesis time. Pre-computed, versioned, pinned to a release. The reads are what the agent does at hypothesis time: single, schema-typed tool calls into pre-computed indices, each returning a structured result with a provenance string. Reasoning happens only at the final step, composing a mechanistic story over the evidence the indices handed it. No inline compute, no retries, no schema drift. Every claim in every report cites the MCP tool that produced it.

This is what we mean by an inference-oriented architecture. The agent's job is to reason over pre-built indices, not to run science. That boundary is what makes the output auditable.

It doesn't fully solve the validity problem on its own. For that we added a prioritization layer: before any hypothesis reaches a human reviewer, we check it against the published literature. Strong published support moves it up the queue. Weak or absent support flags it for closer scrutiny. The scientist sees candidates in rough order of how much the field already believes them, which is the right starting point for deciding whether to commit wet lab resources.

This post is the v0 working example.

---

## What v0 actually does

A typical AMD hypothesis flow looks like this:

1. Start with a GWAS.
2. **Which variants are in high LD?** *(compute: fine-mapping. read: variants, PIPs.)*
3. **Which genes do these variant sets colocalize with?** *(compute: coloc. read: PP4, gene, tissue.)*
4. **Do we see these genes differentially expressed in disease? In what cell types?** *(compute: DE + GSEA. read: log2FC, padj, cell type.)*
5. **What are the functional consequences of the variants?** *(compute: VEP / fold prediction. read: consequence terms, residue, ΔΔG.)*
6. **For which pathways are these genes enriched?** *(compute: enrichment. read: pathway, FDR, overlap.)*
7. **How do all the variants combine to produce the phenotype?** *(compute: none — pure reasoning over everything above.)*
8. **Can we find literature corroboration?** *(compute: literature mining. read: papers that match.)*

A human analyst executes this one step at a time, fixing each break as it happens: cold-start a compute job, wait, retry, paste output into the next step, fix a column mismatch, retry again. When you ask an AI agent to work the same way you inherit every failure mode, process kills, machine timeouts, schema drift, retries, compounded across whatever fanout the workflow has. Six genes × eight steps is forty-eight places to fail.

In v0, every one of those steps is a pre-computed read. You hand the agent an AMD GWAS sumstats file (`samples/amd_fritsche_2016.sumstats.tsv`, GCST003219, Fritsche LG et al. 2016, *Nat Genet*). It reads the `## study_id:` header and calls a registered workflow (`Run_AMD_Hypothesis_v0`). Each step is a single MCP tool call:

| # | Step | MCP server.tool | What's pre-computed |
|---|---|---|---|
| 0 | Read GWAS sumstats | *(local file)* | — |
| 1 | Confirm study in catalog | `jarvis-ot.study_lookup` | Open Targets study table |
| 2 | Fine-mapped credible sets | `jarvis-ot.credible_sets_for_study` | OT credible_set (SuSiE) |
| 3 | Locus-to-gene assignment | `jarvis-ot.l2g_top_genes` | OT L2G predictions (29-feature gradient boosting; coloc evidence baked in as features) |
| 4 | Gene metadata + lead variant | `jarvis-ot.gene_metadata`, `jarvis-ot.lead_variant_for_locus` | OT target + credible_set.locus |
| 5 | **ESM3 fold + variant viz** | `jarvis-esm3.score_target` | Ensembl VEP REST → ESM3 Forge fold → PyMOL render |
| 6 | DE by cell type *(v0 mock)* | `jarvis-indices.query_differential_expression` | Curated AMD scRNA atlas (mock — real backend in v1) |
| 7 | Pathway membership *(v0 mock)* | `jarvis-indices.query_pathway_membership` | Reactome v96 (mock — real backend in v1) |
| 8 | Literature corroboration | `jarvis-paperclip.literature_for_gene` | PaperClip (BM25 + vector over PMC / bioRxiv / medRxiv / arXiv) |
| 9 | Compose mechanistic hypothesis | *(agent reasoning, writes a markdown report)* | — |

Steps 1–8 pass IDs forward with no judgment. Step 9 is the only place reasoning happens. Every claim in the final report cites the MCP tool that produced it.

---

## The six AMD targets

The L2G step returns 20 candidate genes for GCST003219, each with a score and 29-feature SHAP breakdown. The top six, picked by score and biological interpretability, span AMD's two known pillars (complement, lipid) plus the ARMS2/HTRA1 locus.

| Rank | Gene | L2G | UniProt | Lead variant | Consequence | Report |
|---:|---|---:|---|---|---|---|
| 1 | **TNFRSF10A** | 0.966 | O00220 | `8_23225458_G_T` | splice_region_variant | [TNFRSF10A](output/TNFRSF10A_ENSG00000104689.md) |
| 2 | **APOE** | 0.961 | P02649 | `19_44908684_T_C` | **missense C130R (= APOE ε4)** | [APOE](output/APOE_ENSG00000130203.md) |
| 3 | **C9** | 0.961 | P02748 | `5_39327786_G_T` | intron_variant | [C9](output/C9_ENSG00000113600.md) |
| 4 | **CETP** | 0.882 | P11597 | `16_56963437_C_CA` | intron_variant | [CETP](output/CETP_ENSG00000087237.md) |
| 5 | **ARMS2** | 0.863 | P0C7Q5 | `10_122456049_T_C` | intron_variant | [ARMS2](output/ARMS2_ENSG00000254636.md) |
| 6 | **C3** | 0.844 | P01024 | `19_6718376_G_C` | **missense R102G (rs2230199)** | [C3](output/C3_ENSG00000125730.md) |

Two of the six are coding missense variants in genes with decades of literature, **C3 R102G** and **APOE ε4**, both with structural consequences the agent can show, not just describe. The other four are non-coding; the ESM3 step renders the predicted protein and the workflow leans on coloc-as-L2G-feature evidence for gene assignment. The reports are honest about which mode applies to which target.

### What a report looks like

Each per-target report opens with the bottom line: a one-line plain-English hypothesis, then an ASCII flow diagram of the variant-to-AMD chain with the supporting evidence in-line on each connector (SHAP features, Reactome pathway IDs, DE log₂FC and padj values, paper IDs). A "How to verify this evidence" block underneath tells you exactly which MCP tool call produced each citation and how to re-derive it.

The body of the report then layers the detail:

- VEP consequence (with PolyPhen / SIFT predictions for missense)
- L2G SHAP-contributing features, so the reader sees whether distance, coloc, VEP, e2g enhancer-to-gene, or local gene density drove the call
- ESM3 mean pLDDT + pTM, with the variant residue PNG embedded inline
- DE rows for any cell types the v0 atlas covers
- Reactome pathway membership
- Top PaperClip papers with summaries and links
- An agent-composed mechanistic-hypothesis paragraph (written by `claude -p` over the evidence pack at build time; the prose long-form of the ASCII chain at the top)
- Full provenance chain — every claim back to a single MCP tool call

The C3 report shows the R102 side chain in red against the cyan cartoon. The CETP report shows the full β-barrel fold; the variant is intronic so the whole protein renders. The visual mode adapts to what the variant is.

---

## What's real in v0 and what's a stub

The stubs are visible by design. Handwaving them would contradict the whole thesis: if the argument is that inference-oriented architecture requires honest pre-computation, the post can't pretend the mocks are real.

**Real, today:**
- **Open Targets Platform release 2026-03** — `study/`, `credible_set/`, `l2g_prediction/`, `target/` parquets, served via a DuckDB-backed MCP server (`jarvis-ot`). Sizes: 92 MB / 3.8 GB / 530 MB / 81 MB.
- **ESM3 Forge** (`esm3-open-2024-03`) — wrapped in `jarvis-esm3` for fold + InterPro function annotations. Per-protein call: ~13 s for medium proteins, ~25 s for C3 at 1663 aa.
- **Ensembl VEP REST** — variant consequence + PolyPhen/SIFT, cached.
- **UniProt REST** — canonical FASTA, cached.
- **PyMOL open-source headless** — variant residue render to PNG, ~3–10 s.
- **PaperClip** (paperclip.gxl.ai) — BM25 + vector over PMC / bioRxiv / medRxiv / arXiv.

**v0 mocks:**
- **`jarvis-indices.query_differential_expression`** — backed by a small curated AMD scRNA atlas (Orozco LD et al. 2020 *Cell Rep* 30:1246). Good for the top genes; sparse elsewhere. v1 replaces this with a pre-computed DE store over a real single-cell atlas, plus GEO-metadata integration for tissue selection by trait.
- **`jarvis-indices.query_pathway_membership`** — backed by Reactome v96 GMT for the top AMD genes. v1 adds Wikipathways + STRING and binds against the full L2G gene list.

---

## Throughput: human vs agent

A serious analyst working one of these genes from scratch — querying OT, running ESM3, rendering with PyMOL, pulling DE from an atlas, looking up Reactome, surveying literature, writing the synthesis — is in for about five hours per gene if they're fluent, more like a full day if they're learning the tools. Six genes is a focused week of work. The agent does the same six in roughly four minutes warm, seven and a half cold, because every retrieval step is a pre-computed read and the only thing left to do live is reason over the evidence pack:

| Step | Human | Agent (cold) | Agent (warm) |
|---|---:|---:|---:|
| 0. Read GWAS file | 5 min | <0.1 s | <0.1 s |
| 1–3. Study + credible sets + L2G | 30 min | 125 ms | 125 ms |
| 4. Gene meta + lead variant × 6 | 60 min | ~4 s | ~4 s |
| 5a. ESM3 fold × 6 | 90 min | ~3 min | 50 ms |
| 5b. PyMOL variant render × 6 | 30 min | ~30 s | 50 ms |
| 6. DE atlas query × 6 | ~4 h | <0.2 s *(v0 mock)* | <0.2 s *(v0 mock)* |
| 7. Pathway enrichment × 6 | 90 min | <0.2 s *(v0 mock)* | <0.2 s *(v0 mock)* |
| 8. Literature × 6 | ~3.5 h | ~24 s | ~24 s |
| 9. Compose hypothesis × 6 (Claude reasoning) | ~7.5 h | ~3.7 min | ~3.7 min |
| **Total (6 genes)** | **~30 hours** | **~7.5 min** | **~4 min** |
| **Total per gene** | **~5 hours** | **~75 s** | **~42 s** |
| **Speedup per gene vs. human** | — | **~240×** | **~430×** |

Both columns include reads *and* writes. The agent's warm case skips compute via cache — which is the whole point of separating reads from writes. You pay the write cost once, amortize it across every subsequent read. The human re-spends the same time the next time they ask the same question.

At a hundred hypotheses the multiplier diverges further: humans stay at ~5 hours per gene with no amortization, while the agent's per-gene cost approaches retrieval (~5 s) plus composition (~37 s) — and the composition parallelizes trivially across cores. That's the architectural argument restated in throughput terms.

The reasoning step now dominates the agent's runtime. Step 9 — Claude composing the mechanistic-hypothesis paragraph over the evidence pack — eats ~87% of the total. The retrieval substrate, after pre-computation, is essentially free. When the only remaining latency is the model thinking, you're in the right architectural regime.

---

## On MCP latency

Reasonable pushback: if the pitch is "pre-compute aggressively and read fast," doesn't MCP add an IPC layer you don't need? Why not call Python functions directly?

We measured rather than guessed. Median and p95 over 50 iterations of each `jarvis-ot` tool, 10 iterations of each remote call:

| Operation | median | p95 | Source of cost |
|---|---:|---:|---|
| **FastMCP wrapper around a tool call** | **0.96 ms** | 2.04 ms | framework overhead |
| `jarvis-ot.study_lookup` (indexed slim cache) | 0.98 ms | 2.16 ms | hot path |
| `jarvis-ot.credible_sets_for_study` (indexed) | 7.4 ms | 15.1 ms | indexed range scan |
| `jarvis-ot.l2g_top_genes` (materialized join) | 117 ms | 312 ms | indexed scan over 2.8 M rows |
| `jarvis-ot.gene_metadata` (parquet view) | 47.6 ms | 58.6 ms | one call per gene |
| `jarvis-ot.l2g_feature_contributions` | 286 ms | 368 ms | one call per gene |
| `jarvis-ot.lead_variant_for_locus` | 675 ms | 829 ms | one call per gene |
| UniProt REST | 584 ms | 606 ms | network |
| Ensembl VEP REST | 1778 ms | 19844 ms | network + remote, rate-limited tail |
| PaperClip search | 4186 ms | 14522 ms | remote BM25 + vector search |
| PyMOL render PNG (ray-traced) | 3–10 s | — | local CPU |
| ESM3 Forge fold (cold) | 13–30 s | — | remote GPU |

The benchmark script is in the repo (`prototype/scripts/bench_mcp_latency.py`) so the numbers are reproducible.

The FastMCP wrapper costs about a millisecond, statistically indistinguishable from a direct Python call. It's well below the noise floor for anything the agent actually does. A pre-joined DuckDB cache matters more: the first cut of `jarvis-ot` queried raw parquets and `l2g_top_genes` took ~4.7 seconds because every call rescanned 200 + 200 + 10 parquet files. Materializing a slim `(studyId, studyLocusId, chromosome, position, geneId, gene_symbol, l2g_score)` join into a 654 MB indexed DuckDB file dropped the same call to 117 ms, about 40× faster, while keeping struct-heavy columns (the credible-set `locus` variant-membership array, the L2G `features` SHAP struct, the `target.transcripts` blob) on parquet for the calls that need them. Those slower calls run once per gene in the workflow, so 300–700 ms each is fine. The cache build runs in ~90 s on 2 GB of RAM with a memory limit and disk spill (`prototype/scripts/build_ot_cache.py`).

What MCP earns at this scale: process isolation, so ESM3 SDK, DuckDB, PyMOL, and PaperClip can't break each other's startup or runtime; schema-first discovery, so the agent reads tool signatures at session start with no prompt engineering needed; auth boundaries, so `jarvis-esm3` owns the Forge key and `jarvis-paperclip` owns the GXL OAuth token and nothing else touches them; and federation, so when the OT parquet store moves to a memory-rich machine and the ESM3 wrapper moves to a GPU box, the agent doesn't notice.

What MCP doesn't buy: speed in the absolute sense, or batch high-throughput. If your workload is a thousand hypotheses per second on a cluster, skip MCP, ship a single binary with everything in-process, and use gRPC or Arrow Flight at the team boundaries.

The right shape is hybrid. MCP at the agent-facing boundary where discovery, isolation, and auth matter. Direct calls within a service. Spark or Polars on the write side. For v0, one user, one machine, six reports in about four minutes warm, MCP is the right pick and the overhead is invisible. It stays right through multi-user, multi-trait, dozens of concurrent investigations. It stops being right somewhere around millions of hypotheses per second on a cluster, and at that scale there are other problems to solve first.

---

## Why this matters

The test that matters more than throughput numbers: a junior analyst, an MD scientist, and an agent running a thousand investigations in parallel should all execute the same hypothesis flow from the same substrate, without re-running a single compute step. When reads and writes are separated, that becomes possible. The intelligence becomes shared infrastructure. The reasoning becomes the only step that still requires judgment, and because the evidence trail is clean and versioned, you can actually evaluate whether the judgment is good.

That's what makes agent-generated hypotheses trustworthy enough to act on. Not model capability. Infrastructure.

v0 is a demo. v1 will be the same shape, with real DE and real pathways behind the last two stubs, and a `jarvis-mr` Mendelian randomization server alongside.

---

## Try it

Repo at <https://github.com/dorkosaurus/JARVIS_for_bio>. Full workflow at `prototype/mcp_servers/`; AMD demo input at `samples/amd_fritsche_2016.sumstats.tsv`. With Claude Code installed and `.mcp.json` picked up, hand the agent the sumstats file, the workflow ID resolves automatically, and you watch six reports get generated.

Video walkthrough: `[video URL]`.

Pushback welcome, particularly on which compute steps should move next from per-question to pre-computed.

---

*Built on Open Targets, ESM3 (Evolutionary Scale / BioHub), PyMOL, Ensembl, UniProt, PaperClip (GXL), Reactome. The architecture borrows CQRS from the database world. The error of mixing reads and writes is mine alone.*
