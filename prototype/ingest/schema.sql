-- JARVIS prototype: pre-computed indices for AMD hypothesis generation.
-- One SQLite file holds five domain indices + the workflow DAG.
-- Every domain row carries provenance (source dataset + retrieval timestamp)
-- so the agent's final hypothesis can cite where each piece of evidence came from.

PRAGMA foreign_keys = ON;

-- =========================================================================
-- Step 1: credible sets from fine-mapping
-- =========================================================================
CREATE TABLE IF NOT EXISTS credible_sets (
    variant_id        TEXT NOT NULL,           -- rsID, e.g. rs1061170
    gwas_id           TEXT NOT NULL,           -- e.g. GCST90043776 (AMD 2023)
    trait             TEXT NOT NULL,           -- "age-related macular degeneration"
    locus             TEXT NOT NULL,           -- gene-symbol-anchored locus label, e.g. "CFH"
    pip               REAL NOT NULL,           -- posterior inclusion probability
    chrom             TEXT NOT NULL,
    pos               INTEGER NOT NULL,        -- GRCh38
    ref               TEXT,
    alt               TEXT,
    fine_mapper       TEXT NOT NULL,           -- "SuSiE" | "FINEMAP" | "OTG_L2G"
    source_dataset    TEXT NOT NULL,           -- provenance: where row came from
    retrieved_at      TEXT NOT NULL,           -- ISO-8601
    PRIMARY KEY (variant_id, gwas_id)
);
CREATE INDEX IF NOT EXISTS idx_credible_gwas    ON credible_sets(gwas_id);
CREATE INDEX IF NOT EXISTS idx_credible_locus   ON credible_sets(locus);

-- =========================================================================
-- Step 2: colocalization between GWAS variants and molecular QTLs
-- =========================================================================
CREATE TABLE IF NOT EXISTS coloc_results (
    variant_id        TEXT NOT NULL,
    gene_symbol       TEXT NOT NULL,
    gene_ensembl      TEXT,
    tissue_label      TEXT NOT NULL,           -- human-readable, e.g. "retinal pigment epithelium"
    tissue_ontology   TEXT NOT NULL,           -- CL:NNNNN or UBERON:NNNNN
    qtl_type          TEXT NOT NULL,           -- "eQTL" | "pQTL" | "sQTL"
    pp4               REAL NOT NULL,           -- coloc posterior P(shared causal variant)
    qtl_dataset       TEXT NOT NULL,           -- e.g. "GTEx_v8" | "eQTL_Catalogue_r6"
    source_dataset    TEXT NOT NULL,
    retrieved_at      TEXT NOT NULL,
    PRIMARY KEY (variant_id, gene_symbol, tissue_ontology, qtl_type, qtl_dataset)
);
CREATE INDEX IF NOT EXISTS idx_coloc_variant    ON coloc_results(variant_id);
CREATE INDEX IF NOT EXISTS idx_coloc_gene       ON coloc_results(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_coloc_tissue     ON coloc_results(tissue_ontology);

-- =========================================================================
-- Step 3: pathway membership (gene -> pathway)
-- =========================================================================
CREATE TABLE IF NOT EXISTS pathway_membership (
    gene_symbol       TEXT NOT NULL,
    pathway_id        TEXT NOT NULL,           -- e.g. WP:WP4884 or R-HSA-166658
    pathway_name      TEXT NOT NULL,
    pathway_db        TEXT NOT NULL,           -- "WikiPathways" | "Reactome"
    evidence_score    REAL NOT NULL,           -- 0..1; curator-assigned confidence
    n_members         INTEGER,                 -- pathway size (for context-relevance ranking)
    source_dataset    TEXT NOT NULL,
    retrieved_at      TEXT NOT NULL,
    PRIMARY KEY (gene_symbol, pathway_id)
);
CREATE INDEX IF NOT EXISTS idx_pathway_gene     ON pathway_membership(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_pathway_id       ON pathway_membership(pathway_id);

-- Companion table: pathway -> member list (so the agent can ask "what else is in this pathway?")
CREATE TABLE IF NOT EXISTS pathway_members (
    pathway_id        TEXT NOT NULL,
    gene_symbol       TEXT NOT NULL,
    PRIMARY KEY (pathway_id, gene_symbol)
);
CREATE INDEX IF NOT EXISTS idx_pathway_members_pid ON pathway_members(pathway_id);

-- =========================================================================
-- Step 4: gene programs from single-cell decomposition (NMF / topic modeling)
-- =========================================================================
CREATE TABLE IF NOT EXISTS gene_programs (
    program_id          TEXT NOT NULL,         -- e.g. "RETINA_RPE_NMF_07"
    program_label       TEXT NOT NULL,         -- human-readable theme, e.g. "complement activation"
    cell_type_label     TEXT NOT NULL,
    cell_type_ontology  TEXT NOT NULL,         -- CL:NNNNN
    atlas_dataset       TEXT NOT NULL,         -- e.g. "Voigt_2019_choroid_scRNA"
    method              TEXT NOT NULL,         -- "NMF" | "cNMF" | "topic_model"
    coherence           REAL,                  -- program coherence score
    source_dataset      TEXT NOT NULL,
    retrieved_at        TEXT NOT NULL,
    PRIMARY KEY (program_id)
);

-- gene -> program membership (top-loaded genes per program)
CREATE TABLE IF NOT EXISTS gene_program_members (
    program_id        TEXT NOT NULL,
    gene_symbol       TEXT NOT NULL,
    loading           REAL NOT NULL,           -- NMF loading or topic weight
    rank              INTEGER NOT NULL,        -- 1 = top-loaded
    PRIMARY KEY (program_id, gene_symbol),
    FOREIGN KEY (program_id) REFERENCES gene_programs(program_id)
);
CREATE INDEX IF NOT EXISTS idx_gpm_gene  ON gene_program_members(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_gpm_prog  ON gene_program_members(program_id);

-- =========================================================================
-- Step 5: differential expression (case vs control, per cell type)
-- =========================================================================
CREATE TABLE IF NOT EXISTS differential_expression (
    gene_symbol         TEXT NOT NULL,
    cell_type_label     TEXT NOT NULL,
    cell_type_ontology  TEXT NOT NULL,
    comparison          TEXT NOT NULL,         -- e.g. "AMD_vs_control"
    log2fc              REAL NOT NULL,
    padj                REAL NOT NULL,
    method              TEXT NOT NULL,         -- "DESeq2" | "edgeR" | "MAST"
    atlas_dataset       TEXT NOT NULL,
    n_case              INTEGER,
    n_control           INTEGER,
    source_dataset      TEXT NOT NULL,
    retrieved_at        TEXT NOT NULL,
    PRIMARY KEY (gene_symbol, cell_type_ontology, comparison, atlas_dataset)
);
CREATE INDEX IF NOT EXISTS idx_de_gene    ON differential_expression(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_de_celltype ON differential_expression(cell_type_ontology);

-- =========================================================================
-- Workflow DAG: encodes the canonical hypothesis-generation flow
-- =========================================================================
CREATE TABLE IF NOT EXISTS workflows (
    workflow_id    TEXT PRIMARY KEY,
    objective      TEXT NOT NULL,              -- natural-language objective this flow serves
    description    TEXT
);

CREATE TABLE IF NOT EXISTS workflow_steps (
    workflow_id    TEXT NOT NULL,
    step_id        TEXT NOT NULL,              -- e.g. "step_1_credible_sets"
    step_order     INTEGER NOT NULL,
    step_name      TEXT NOT NULL,
    description    TEXT NOT NULL,              -- what the agent should do at this step
    tool_name      TEXT NOT NULL,              -- MCP tool to call
    input_schema   TEXT NOT NULL,              -- JSON: tool input shape
    output_schema  TEXT NOT NULL,              -- JSON: expected output shape
    depends_on     TEXT,                       -- comma-separated step_ids
    PRIMARY KEY (workflow_id, step_id),
    FOREIGN KEY (workflow_id) REFERENCES workflows(workflow_id)
);
