# Drug discovery has a reasoning problem (building JARVIS v1)

Biological data has scaled past any single human's capacity to reason over it.  Frontier models allow us to break through human reasoning limitations but the infrastructure we have was not designed for agents.

We need an "inference-oriented architecture" that supports agentic-driven workflows at scale.  The key change we need to make:  separate reads from writes. Pre-compute all the data we need to "flow" over questions (e.g. fine-mapping, colocalization, fold prediction, enrichment, and literature validation). 

With such an architecture, I was able to form six AMD mechanistic hypothesis reports along with the evidence that underlies them (C3 R102G, APOE ε4, C9, CETP, ARMS2, TNFRSF10A) in under four minutes.  Generating such reports would take many, many hours for humans to create.  And trying to leverage agents to do this without the pre-computed data would lead to frustration as agents fail because they are marrying compute + inference.   

Write-up + reports + reproducible benchmark: [GITHUB-README-URL]

90-second demo: [VIDEO-URL]

Which compute step should move next from per-question to pre-computed?

#drugdiscovery #AIxBiology #targetdiscovery #BuildingJARVIS

