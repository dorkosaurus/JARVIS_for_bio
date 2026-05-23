# Enabling reason at scale
### An architecture to enable frontier models to reason over all our data

---

## The token economics problem

Biological data is large. Across genomics, proteomics, high-content screening, and the growing diversity of organisms we characterize, biological data generation is at least exabyte-scale.

After dimensionality reduction (transformations from raw instrument output to the data we actually reason on) we can reduce this "raw" data by orders of magnitude. Exabyte-scale raw data will process down to petabyte-scale data we can reason over.

But no human can reasonably reason over petabyte-scale data. So we need machines, and increasingly, we want to use AI.

The naive approach (which many have tried) is to point AI at biological data and say "Go". But that approach is bound to fail because it hits the token economics barrier.

A token is approximately 4 bytes. A single petabyte of data requires 250 trillion tokens to process. At current pricing, passing that data through a frontier model will cost between $12 million and $1.25 billion, depending on the model.

This is a structural barrier we are hitting — we can't optimize our way out of this. Prices would have to fall by orders of magnitude before naive tokenization becomes viable.

---

## The problem with current architectures

The key challenge we face are the architectures we use to store and use data in biology. They are too _data-centric_ and not _inference-centric_.

Let's use a simple example: the process by which a scientist reasons over data to form a therapeutic hypothesis in the context of AMD (Age-related Macular Degeneration):

1. She starts with a GWAS and identifies credible variant associations.
2. She finds a trans-colocalized pQTL linking variant A to gene C, a member of the Wnt signaling pathway.
3. She finds single-cell evidence in retinal pigment epithelial cells showing that gene C is part of the gene program regulating Wnt signaling.
4. She forms the hypothesis that AMD in RPE cells may be driven by dysregulation of the Wnt signaling pathway.

Each step in the investigator's reasoning chain requires a separate computational act:

- Credible set identification requires running fine-mapping software against GWAS summary statistics.
- Colocalization requires executing a statistical model across two datasets simultaneously.
- Cell type expression requires retrieving count matrices, normalizing them, and running differential expression analysis.

Each of these operations is expensive, bespoke, and performed on demand: every time, by every investigator, independently.

This is a data-centric architecture in practice. Data is stored in its rawest useful form and computation happens at query time.

It made sense when investigators were few, datasets were small, and queries were infrequent. It doesn't make sense when we want AI to reason across the entire evidential space of biology continuously and at scale. We'd need AI to _assemble_ all of this data in real time and would quickly hit the token economics barrier.

---

## An inference-oriented architecture

If we want AI to reason over biological data at scale, we need a different architecture that maximally assembles and compresses this data so that the economics become more sensible.



**Components:**

- **A graph database contains a taxonomy of common analytic workflows.** The example provided above gives a canonical human genetics workflow. There are many others. Each node in a DAG corresponds to a task that a scientist performs in pursuit of an objective.

- **Real-time ontology assignment from queries.** Pre-computed data indices should be linkable via ontologies but queries might not be. We need a real-time mechanism to take a query string and transform string components into ontology terms that enable robust retrieval from indices. For example, the query might provide "RPE" and the ontology expander should transform that into "Retinal Pigment Epithelial Cell" (CL:00025860). The Ontology Access Kit (OAK) from the Monarch Initiative provides a nice example of how to do this.

- **Biological foundation models.** Foundational models massively integrate and compress biological data into formats that both humans and AI can reason over. Models like ESM, Evo, and Geneformer internalize biological structure from their data during pretraining. They already "know" what a splice site or a protein domain looks like. You don't tokenize the evidence, you tokenize the question.

- **Pre-computed indices for all the data used by different tasks.** Each node in an analytic workflow will query for one or more data types. For this to scale, we need to precompute all the data we need for query efficiency. This is a clear _buy_ opportunity. We need providers like the team at Devano AI who are working on assembling such indices that companies can buy and add their internal data to.

- **Embedding indexes over biological meaning.** Keywords retrieve what was labeled. Embeddings retrieve what is relevant. For each modality (cell states, protein sequences, perturbation responses, chemical structures) we need dense vector indexes built over biological meaning rather than metadata. A query for RPE cells under Wnt dysregulation surfaces relevant biological states across the entire atlas, including datasets that never explicitly studied AMD.

Finally, not in the graph but a prerequisite for production work: **we need a database that contains per-task benchmarks and evaluations.** Knowing the workflows isn't enough. We need to formally define rubrics that define "good enough" for each task so we can assess the quality of AI decision making and enable realistic financial planning for AI utilization.

---

Will this work? I'll find out in the next two weeks — I'm building a prototype to instantiate this architecture.
