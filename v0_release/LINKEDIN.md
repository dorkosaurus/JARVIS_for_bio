# Building JARVIS: an inference-oriented architecture to enable agent-scale reasoning

Some fields have more questions than humans can ever ask and answer.

* In biology: gene × cell × variant × tissue × pathway.
* In materials science: composition × structure × processing × property.

The hypotheses formed from these questions scale combinatorially. 

AMD alone: 20 loci × 6 candidate genes × 10 protein states × 5 pathways × 4 cell types ≈ 24,000 combinations. And that's one disease.  The GWAS catalog has millions of variants across thousands of diseases.

Now AI agents can finally explore that space. But the infrastructure underneath wasn't built for them. It was built for humans asking one question at a time. Asking an agent to use this infrastructure sets the agent up for failure.

The fix: an inference-oriented architecture. Separate reads from writes. Pre-compute the heavy lifting (fine-mapping, colocalization, fold prediction, enrichment, literature) ahead of time. Minimize reads and where you need computation, make sure it's low latency.   

Such an architecture sets up agents for success and also enables provenance and reproducibility since the agent is reasoning over versioned, read-only data.  

I shipped a v0 release of this architecture and used it to generate six agentically-created AMD mechanistic hypotheses (C3 R102G, APOE ε4, C9, CETP, ARMS2, TNFRSF10A) in under four minutes using low-end Digital Ocean machines. The same work would likely take a human team many, many hours.

Write-up + reports + reproducible benchmark in the comments section. 

I welcome feedback and discussion!  Also, which compute step in your field should move next from per-question to pre-computed?

#AIxScience #agentinfrastructure #drugdiscovery #BuildingJARVIS
