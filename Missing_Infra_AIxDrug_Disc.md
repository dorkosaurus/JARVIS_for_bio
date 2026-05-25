# The missing infrastructure for AI-driven drug discovery

## The role that mechanistic hypotheses play in drug discovery

Scientists in drug discovery spend a lot of their time forming and validating mechanistic hypotheses.  Mechanistic hypotheses describe how molecular perturbations exert their effects across levels of biological organization to produce observed phenotypes. 

There are three primary mechanistic hypotheses in drug development: 

* Therapeutic hypotheses:  these are hypotheses that describe how biological perturbations lead to observed disease states.
* Tox hypotheses: these are hypotheses on how perturbations introduced lead to hazards observed by toxicologists.  
* Mechanism of action:  these are hypotheses that describe how a therapeutic specifically modulates a biology to produce the intended or observed effect.

## How scientists form mechanistic hypotheses

Scientists form mechanistic hypotheses by starting with a "finding" and integrating it across different kinds of data.  When the data line up, it starts telling them a story about the underlying biology.   

To make this concrete, here's a simple example: the process by which a scientist reasons over data to form a therapeutic hypothesis in the context of AMD (Age-related Macular Degeneration):

1. She starts with a GWAS and identifies credible sets for variant associations.
2. She finds a trans-colocalized pQTL linking variant A to increased expression of gene C, a member of the Wnt signaling pathway.
3. She finds single-cell evidence in retinal pigment epithelial cells showing that gene C is part of the gene program regulating Wnt signaling.
4.  She finds differential expression evidence that suggests that these gene programs are uniformly upregulating in AMD patients vs controls. 
5. She then interprets that AMD in RPE cells may be driven by dysregulation of the Wnt signaling pathway by upregulating gene C. 
6. This leads to a therapeutic hypothesis:  downregulation of gene C might ameliorate AMD symptoms via better regulating the WNT signalling pathway in RPE cells.

## Mechanistic hypotheses are analytic flows

Mechanistic hypotheses are just analytic flows and they follow patterns.  In the example above, the scientist used a canonical GWAS analysis pattern:

1. Create credible sets of variants within the GWAS.
2. Perform causal gene analysis (including pQTL colocalization).
3. Identify pathways it plays a role in — these become gene sets we want to interrogate.
4. Identify gene programs from single-cell expression data where we can find these gene programs.
5. Perform differential expression between case/controls to determine if pathways are differentially expressed between people who have AMD and people who do not.
6. Form the hypothesis.

## How AI can scale our ability to perform analytic workflows

It's a lot of work for scientists to go through all of these data.  If we could find a way for AI to reliably execute all of these analytic flows at scale (interrogating every finding across all relevant and available data), we could explore many such hypotheses and give them to scientists for validation.  

Note, this would create a new problem to solve - hypothesis validation and prioritization but this is a good problem to have. 

## Why AI can't do this for us today

The key problem is that agents don't have the environment they need to successfully perform these tasks.  Specifically:

1.  They lack well described analytic flows that are programmatically accessible and executable. 
2.  They lack infrastructure designed to help them efficiently reason over this data as they execute these analytic flows.  Since biological evidence is at least petabyte scale, this could burn **trillions** of tokens at very high cost.  

To make this concrete, here's an overview on reasons why agents might fail to run the workflow described above:

**Step 1: Credible set generation.** An agent handed a raw GWAS summary statistic file needs to run fine-mapping software (SuSiE, FINEMAP) to generate credible sets. These tools require compute infrastructure, reference LD panels matched to the study population, and parameter choices that affect downstream results. No agent can invoke these reliably without a pre-configured execution environment and validated reference data. The agent doesn't fail because it doesn't know what to do — it fails because the tools aren't callable.

**Step 2: Causal gene assignment and pQTL identification.** Colocalization against a pQTL dataset requires knowing which pQTL dataset to use, which tissues are relevant, and what statistical thresholds are appropriate. The eQTL Catalogue contains hundreds of datasets. An agent querying it naively will either retrieve too much (burning tokens) or make arbitrary choices about which dataset is relevant without the context to justify them.

**Step 3: Pathway identification.** Mapping a gene to pathways is tractable — databases like WikiPathways and Reactome are queryable. But the output is a gene set, and the agent now needs to decide which gene sets are biologically meaningful in the context of AMD specifically, not just which ones contain gene ABC. That judgment requires context the agent doesn't have unless it's been pre-loaded.

**Step 4: Gene program identification from single-cell data.** This is where the token economics becomes insurmountable without pre-computed infrastructure. Identifying gene programs from a single-cell atlas requires access to count matrices, dimensionality reduction, and gene program analysis (NMF, topic modeling). The relevant atlas for AMD — cells from retinal tissue across AMD and control donors — may be hundreds of gigabytes. An agent cannot retrieve and reason over this at query time.

**Step 5: Differential expression.** Running differential expression between AMD and control RPE cells requires the same count matrices, a statistical framework (DESeq2, edgeR), and choices about covariates, normalization, and multiple testing correction. Not callable without pre-configured infrastructure and pre-selected datasets.

**Step 6: Hypothesis formation.** This is the one step an agent can do well — if it has been given the outputs of steps 1–5 in a compact, well-structured form. The hypothesis formation step is a reasoning task. Every step before it is a data engineering task. Agents are being asked to do both, and the data engineering is what breaks them.

## An analysis-centric architecture

The key architectural insight is that every step in a biological analytic workflow is a data retrieval problem, not a compute problem. The compute happened — or should have happened — at ingest time. The agent's job is to retrieve pre-computed results, reason over them, and decide what to retrieve next.

This requires three things:

1. A programmatically accessible taxonomy of analytic flows. The AMD workflow above is one canonical pattern. There are many others — starting from expression data, from compound screens, from rare disease phenotypes, from protein structure. Each flow needs to be encoded as a DAG where every node has a well-defined input schema, output schema, and a pointer to the infrastructure that serves it. The agent follows the flow. It doesn't invent it.
2. Pre-computed indices for every data type used by every step. Each arrow in the analytic flow represents a query. That query should return a pre-computed answer in milliseconds, not trigger a compute job that runs for hours. Credible sets, colocalization results, pathway memberships, gene programs, differential expression statistics — all of it materialized at ingest, versioned, and queryable via a lightweight API. The agent never sees raw data. It sees structured, provenance-tagged answers.
3. Ontology-grounded retrieval. Biological queries are almost always under-specified. "Retinal pigment epithelial cells" needs to resolve to every dataset containing that cell type regardless of how the original authors labeled it, before a single record is retrieved. This requires real-time ontology expansion — transforming query strings into controlled vocabulary terms — before hitting any index.

## Prototyping an analysis-centric architecture

I wanted to create a v0 prototype of this architecture to demonstrate the pieces and user experience.  [This document](./missing_infra_prototype/README.md) describes it in detail.
