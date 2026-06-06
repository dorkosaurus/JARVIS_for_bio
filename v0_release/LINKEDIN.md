# Building JARVIS: an inference-oriented architecture to enable agent-scale reasoning

Some fields have more questions than humans can ever ask and answer.

* In biology: gene × cell × variant × tissue × pathway.
* In materials science: composition × structure × processing × property.

The hypotheses scale combinatorially. AMD alone: 20 loci × 6 candidate genes × 10 protein states × 5 pathways × 4 cell types ≈ 24,000 combinations. And that's one disease — the GWAS catalog has millions of variants across thousands.

Now AI agents can finally explore that space. But the infrastructure underneath wasn't built for them. It was built for humans asking one question at a time. Asking an agent to use it sets the agent up for failure.

The fix: an _inference-oriented architecture_. Separate reads from writes. Pre-compute the heavy lifting (fine-mapping, colocalization, fold prediction, enrichment, literature) ahead of time. The agent only reads and reasons. Every number in the report traces back to a pre-computed source — forensics by default.

With this architecture, I generated six AMD mechanistic hypothesis reports — C3 R102G, APOE ε4, C9, CETP, ARMS2, TNFRSF10A — in under four minutes using low-end Digital Ocean machines. The same work takes a human team many, many hours.

Write-up + reports + reproducible benchmark: https://github.com/dorkosaurus/JARVIS_for_bio/blob/main/v0_release/README.md

Which compute step in your field should move next from per-question to pre-computed?

#AIxScience #agentinfrastructure #drugdiscovery #BuildingJARVIS
