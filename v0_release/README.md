# An inference-oriented architecture for agent-scale reasoning

Some fields have more questions than humans can ever ask and answer.

* In biology: gene × cell × variant × tissue × pathway.
* In materials science: composition × structure × processing × property.

Anywhere the question space is built from combinations of things, the hypotheses that *could* be formed scale combinatorially. As a result, most questions don't get asked or answered and most hypotheses don't get formed.

What's changed is that AI agents can finally explore hypothesis space at the scale the combinatorics demand.

But agents are running on top of infrastructure designed for humans asking one question at a time: fine-mapping a GWAS locus, folding a protein, running pathway enrichment, mining the literature.

At each step, humans perform both compute and read operations. A geneticist takes a GWAS and runs a fine-mapping method to find credible sets. The method is slow. It times out. It gets killed. The human re-runs it and moves on.

Agents can't work this way. And we are asking them to.  We are setting them up for failure.

The arithmetic is unforgiving. A serious AMD target-discovery effort might touch twenty GWAS loci. Each one needs fine-mapping, colocalization, differential expression, variant consequence annotation, pathway enrichment, and literature synthesis before you can say anything mechanistic about it.

And the space to explore isn't twenty. It's twenty loci × six candidate genes per locus × ten protein states × five pathways × four AMD-relevant cell types. Around 24,000 combinations. Most never get explored.

And that's just AMD. The GWAS catalog has millions of variants across thousands of diseases. Each one runs the same combinatorial fan-out.

Throughput isn't the only problem. Agents lack _scientific_ judgement — they don't have the decades of training or the nuance an average scientist has. They don't know when a colocalization PP4 of 0.6 is compelling versus marginal, or that ARMS2 has been contested in the AMD literature for fifteen years. They generate confidently regardless. We're shifting the bottleneck from hypothesis generation to hypothesis validation.

There's a third problem: forensics. The agent does both the compute and the reasoning. When a number lands in the report, you can't tell whether the tool produced it or the agent inferred it. An MD scientist asking "where did this PP4 of 0.6 come from?" gets no clean answer. You can't fix what you can't trace.

Throughput and forensics point to the same architectural fix: separate reads from writes. Pre-computation enables both.

Writes happen once, in batch, ahead of hypothesis time. Fine-mapping, colocalization, fold prediction, enrichment. Pre-computed. Versioned. Pinned to a release.

Reads are what the agent does at hypothesis time. Single, schema-typed tool calls into pre-computed indices. Each one returns a structured result with a provenance string.

Reasoning happens only at the final step. The agent composes a mechanistic story over the evidence the indices handed it. No inline compute. No retries. No schema drift. Every claim in every report cites the MCP tool that produced it.

We need a name for an environment that lets agents reason at scale. One doesn't exist. So we're calling it an _inference-oriented architecture_. The agent reasons over pre-built indices. It doesn't perform significant computation, because the architecture is built for low-latency operations. Every number in the report now traces back to a single pre-computed source (enabling forensics).

Validation needs a different fix. Before any hypothesis reaches a scientist, we check it against the published literature. Strong support moves it up the queue. Weak or absent support flags it for review. The scientist sees candidates ranked by what the field already supports.

I've built a v0 version of an inference-oriented architecture — described below.

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

What MCP earns at this scale:

* **Process isolation.** ESM3 SDK, DuckDB, PyMOL, and PaperClip can't break each other's startup or runtime.
* **Schema-first discovery.** The agent reads tool signatures at session start. No prompt engineering needed.
* **Auth boundaries.** `jarvis-esm3` owns the Forge key. `jarvis-paperclip` owns the GXL OAuth token. Nothing else touches them.
* **Federation.** When the OT parquet store moves to a memory-rich machine and the ESM3 wrapper moves to a GPU box, the agent doesn't notice.

What MCP doesn't buy: speed in the absolute sense, or batch high-throughput. If your workload is a thousand hypotheses per second on a cluster, skip MCP, ship a single binary with everything in-process, and use gRPC or Arrow Flight at the team boundaries.

The right shape is hybrid. MCP at the agent-facing boundary where discovery, isolation, and auth matter. Direct calls within a service. Spark or Polars on the write side. For v0, one user, one machine, six reports in about four minutes warm, MCP is the right pick and the overhead is invisible. It stays right through multi-user, multi-trait, dozens of concurrent investigations. It stops being right somewhere around millions of hypotheses per second on a cluster, and at that scale there are other problems to solve first.

---

## Why this matters

The test that matters more than throughput numbers: a junior analyst, an MD scientist, and an agent running a thousand investigations in parallel should all execute the same hypothesis flow from the same substrate. No-one re-runs compute.

When reads and writes are separated, that becomes possible. The intelligence becomes shared infrastructure. Reasoning is the only step that still requires judgment. And because the evidence trail is clean and versioned, you can actually evaluate whether the judgment is good.

That's what makes agent-generated hypotheses trustworthy enough to act on. Not model capability. Infrastructure.

The shape of the fix isn't unique to biology. Anywhere the question space is combinatorial — materials, chemistry, climate, code — the same separation of reads from writes is what turns an agent from a frustrated single-question worker into something that can actually canvas the space.

v0 is a demo. v1 will be the same shape, with real DE and real pathways behind the last two stubs, and a `jarvis-mr` Mendelian randomization server alongside.

---

## Try it

Repo at <https://github.com/dorkosaurus/JARVIS_for_bio>. Full workflow at `prototype/mcp_servers/`; AMD demo input at `samples/amd_fritsche_2016.sumstats.tsv`. With Claude Code installed and `.mcp.json` picked up, hand the agent the sumstats file, the workflow ID resolves automatically, and you watch six reports get generated.

Video walkthrough: `[video URL]`.

Pushback welcome, particularly on which compute steps should move next from per-question to pre-computed.

---

*Built on Open Targets, ESM3 (Evolutionary Scale / BioHub), PyMOL, Ensembl, UniProt, PaperClip (GXL), Reactome. The architecture borrows CQRS from the database world. The error of mixing reads and writes is mine alone.*
