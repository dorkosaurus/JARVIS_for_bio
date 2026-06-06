# Social posts

Copy-paste ready. Substitute `[substack-url]` and `[video-url]` once posted.

---

## LinkedIn (5 lines)

> Drug discovery needs to separate reads from writes.
>
> Mixing compute and lookup inside every analytic step is what makes agent-driven hypothesis generation fall over. Pre-compute the reads, ship them as MCP tools, and let the agent reason only at the last step.
>
> I built a v0 that does this end-to-end for an AMD GWAS — six mechanistic hypothesis reports in ~60 seconds, generated from `samples/amd_fritsche_2016.sumstats.tsv` through Open Targets L2G, ESM3 fold + variant viz, and PaperClip literature. Every claim cites the MCP tool that produced it.
>
> The full write-up is here: [substack-url]. Demo video here: [video-url].
>
> Pushback welcome — particularly on which compute steps should move from per-question to pre-computed next.

---

## LinkedIn (alt — shorter, sharper opener)

> AI agents can't scale biological hypothesis generation while every analytic step tangles compute with lookup.
>
> CQRS solved this for databases two decades ago: separate the reads from the writes. The reads (fine-mapping, coloc, fold prediction, enrichment) get pre-computed and served as MCP tools. The agent reads them at hypothesis time and reasons only at the end. No timeouts, no schema drift, full provenance.
>
> v0 runs end-to-end for an AMD GWAS: six mechanistic per-target reports in ~60s. Hero example is the C3-R102G missense — Open Targets L2G picks the gene, ESM3 folds it, PyMOL renders the variant residue, PaperClip surfaces the corroborating literature.
>
> Full write-up: [substack-url] · Video walkthrough: [video-url]
>
> Built on Open Targets, ESM3 / BioHub, Ensembl VEP, UniProt, PaperClip (GXL), Reactome.

---

## X / Twitter (1 line, ~280 chars)

Option A — thesis-first:

> Drug discovery needs to separate reads from writes. v0 of an inference-oriented architecture for agent-driven hypothesis generation — AMD GWAS to 6 mechanistic reports in ~60s, end-to-end on MCP. Write-up: [substack-url] · Demo: [video-url]

Option B — demo-first:

> Hand a Claude agent an AMD GWAS, get 6 mechanistic hypothesis reports back — Open Targets L2G + ESM3 fold + variant viz + PaperClip literature, all over MCP, every claim citing the tool that produced it. v0 of an inference-oriented bio stack: [substack-url] · Video: [video-url]

Option C — concrete-finding-first:

> Built a v0 inference-oriented stack for agent-driven drug discovery. AMD GWAS in → 6 mechanistic per-target reports out (incl. C3-R102G with ESM3 fold + PyMOL variant render), end-to-end on MCP servers in ~60s. Read/write separation, CQRS style. [substack-url] · [video-url]

Pick whichever frames the audience you want.
