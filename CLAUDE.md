# JARVIS-for-bio prototype

Inference-oriented prototype: Claude Code orchestrates analyses by calling
MCP-exposed indices and following pre-registered analytic workflows.
Design doc: `Missing_Infra_AIxDrug_Disc.md`. Implementation notes:
`missing_infra_prototype/README.md`.

## MCP servers (auto-loaded via .mcp.json)

- **jarvis-workflow** — `list_workflows`, `get_workflow`,
  `find_workflow_for_objective`. Call this **first** when given an
  analytical objective. Workflows are pre-registered analytic DAGs; follow
  them rather than improvising.
- **jarvis-indices** — 8 tools over the 5 indices (credible sets, coloc,
  pathways, gene programs, differential expression). Every tool returns a
  `provenance` field; cite it in your final answer.

## Sample GWAS datasets

`missing_infra_prototype/samples/*.sumstats.tsv` are GWAS summary statistics files
(top-loci preview format, GWAS Catalog harmonized style). Each file has
`##` comment headers at the top — the line `## study_id: <id>` is the
key that maps to a study already indexed in `credible_sets`. When the
user references one of these files or says "this GWAS":

1. Read the file. Extract `study_id` from the `## study_id:` header.
2. Use that ID to query the jarvis-indices server (fine-mapping, coloc,
   pathway, programs, DE are all **pre-computed** — do not re-run them).
3. The TSV body shows ~10 top loci for context; do not claim to have
   parsed the full summary statistics or re-fine-mapped from them.

## Conventions

**1. Workflow before action.** When the user gives an analytical objective
(e.g., "what are the implications of this GWAS?"), the first move is:

```
jarvis-workflow.find_workflow_for_objective(objective_substring="GWAS")
jarvis-workflow.get_workflow(workflow_id=<id from above>)
```

**2. Plan before execution.** After retrieving the workflow, present the
step list to the user as a short numbered plan (1 line per step: the step
name + which tool it will call + what it operates on). Then ask:

> Shall I proceed?

Wait for explicit approval ("yes", "go", "proceed") before calling any
`jarvis-indices` tool. If the user says no or asks to change something,
revise the plan and re-confirm.

**3. Provenance over invention.** Cite the `provenance` string returned by
each index tool in your final hypothesis. If an index returns no rows for
a query, say so honestly — do not fabricate biology that isn't backed by a
tool result. Surface ambiguity (e.g., causal-gene uncertainty at a locus)
rather than picking by fiat.

**4. Reasoning is the final step.** The first 5–6 steps of a workflow are
retrieval. Reasoning happens only in the final step, over the compact
evidence the prior retrievals returned. Don't reason ahead of the data.
