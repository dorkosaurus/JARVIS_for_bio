# JARVIS for Biology

> *Helping biology move faster by expanding hypothesis space, accelerating experimentation, and increasing the rate of information flow across silos.*

---

## The problem

Biology is generating data faster than we can reason over it. Across genomics, proteomics, high-content screening, and the growing diversity of organisms we characterize, yearly data generation is at least exabyte-scale. After dimensionality reduction, that's multiple petabytes of evidence accumulating annually across modalities, organisms, and biological contexts.

No human can reason over petabyte-scale evidence. And the naive approach — pointing a frontier AI model at biological data and saying "Go" — hits a structural barrier immediately. A single petabyte of data requires 250 trillion tokens to process. At current pricing, that costs between $12 million and $1.25 billion for a single inference pass. This isn't a cost optimization problem. It's an architectural one.

At the same time, biology is siloed. The geneticist studying AMD doesn't know about the structural biologist who characterized the relevant protein. The chemist who screened against the target doesn't know about the single-cell atlas that shows its expression in the relevant cell type. Information exists. It just doesn't flow.

The result: biology moves slower than it should. Hypotheses that could be formed in days take months. Experiments that could be avoided are run. Connections that could accelerate programs go unmade.

---

## The vision

JARVIS for biology is an AI system that operates as a true scientific collaborator — one that can explore hypothesis space at scale, help design and run experiments, and connect evidence across the silos that fragment biological knowledge.

It is organized around three capabilities:

### 1. Explore
*Expand the hypothesis space available to investigators.*

Most hypotheses never get formed because no single investigator can hold all the relevant evidence in mind simultaneously. JARVIS explores hypothesis space systematically — traversing canonical analytic workflows across genetics, expression, chemistry, phenotype, and structure — to surface connections that would otherwise remain invisible.

This requires an inference-oriented data architecture: pre-computed indices, ontology-grounded retrieval, biological foundation models, and embedding indexes over biological meaning. The goal is to make petabyte-scale evidence queryable at the token economics of a conversation. See [Enabling Reason at Scale](./Enabling_reason_at_scale.md) for the full architectural argument.

### 2. Experiment
*Accelerate the transition from hypothesis to evidence.*

A hypothesis is only valuable if it can be tested. JARVIS connects hypothesis generation to experimental execution — helping investigators design experiments, select appropriate model systems, anticipate confounds, and close the loop between prediction and observation. In closed-loop settings, it drives autonomous experimentation cycles without waiting for human intervention at every step.

### 3. Connect
*Increase the rate of information flow across biological silos.*

Biology's knowledge is fragmented across labs, institutions, modalities, and time. JARVIS acts as connective tissue — surfacing relevant work from adjacent fields, linking findings across experimental contexts, and ensuring that evidence generated in one corner of biology becomes available to investigators working in another.  For example, a program team might be working on creating a therapeutic which upregulates expression of a gene because downregulation is thought to create disease state.  In a different part of the org, a human genetics scientist might find that downregulation does __not__ lead to disease state. How can we rapidly share this finding with the program team?  this requires real time contextual monitoring, understanding, and notification.

---

## Architecture

The inference-oriented data layer that enables hypothesis exploration at scale is described in detail in [Enabling Reason at Scale](./Enabling_reason_at_scale.md). Its core components are:

| Component | Purpose |
|---|---|
| Analytic workflow graph | Taxonomy of canonical investigative flows encoded as a DAG |
| Real-time ontology assignment | Transforms query strings into controlled vocabulary terms before retrieval |
| Biological foundation models | Compress biological evidence into queryable representations (ESM, Evo, Geneformer) |
| Pre-computed indices | Materializes task-specific data at ingest rather than query time |
| Embedding indexes | Dense vector retrieval over biological meaning across all modalities |
| Per-task benchmarks | Quality rubrics and token budgets for each step in each analytic workflow |

---

## Current status

This project is under active development. Current work includes:

- **Inference-oriented architecture prototype** — instantiating the components described in [Enabling Reason at Scale](./Enabling_reason_at_scale.md) against a canonical human genetics workflow (GWAS → colocalization → cell type expression → hypothesis)
- **Closed-loop experimentation** — autonomous experimental cycles with Opentrons robotics and Bayesian optimization, active in partnership with Monomer Bio and JBEI
- **Surveillance architecture** — contextual notification system that monitors evidence sources for changes relevant to active hypotheses. Where the inference layer answers questions you knew to ask, surveillance answers questions you didn't know to ask yet — watching public repositories, preprint servers, and data releases and routing signal to investigators based on the context of their current work, not just keyword matches

---

## Why this matters

The rate of biological discovery is a function of how fast investigators can move through hypothesis space, how quickly hypotheses can be tested, and how efficiently evidence flows across the people generating it. JARVIS is designed to accelerate all three.

The goal is not to replace the biologist. It is to give her capabilities she has never had before — the ability to reason over all the evidence, not just the fraction she has time to read; to test hypotheses faster than any single lab could; and to benefit from discoveries made in adjacent fields the moment they become relevant to her work.

Biology should move faster. This is how.

---

## Getting started

*Documentation in progress. Check back as the prototype develops.*

---

## Contact

Vivek Ramaswamy | [Evolution Engines](https://evolutionengines.ai)
