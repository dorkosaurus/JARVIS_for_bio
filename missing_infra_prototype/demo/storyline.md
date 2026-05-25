# AMD demo storyline

This document specifies the biological story the v0 prototype is built around — the
locus, the genes, the cell type, and the hypothesis the agent should be able to reach
after traversing the workflow. Every row in the indices traces back to one of these.

## Anchor disease

**Age-related macular degeneration (AMD)** — MONDO:0008056. Leading cause of
irreversible vision loss in adults over 50 in developed countries.

## Anchor GWAS

**Fritsche et al. 2016** (Nat Genet 48:134) is the canonical large AMD GWAS:
- 16,144 advanced AMD cases vs 17,832 controls
- 52 independent variants at 34 loci reach genome-wide significance
- Open Targets Genetics study ID: **GCST003473**
- (Newer studies exist but this one has the cleanest public credible sets and is
  the reference everyone cites.)

## Anchor loci (primary)

| Locus  | Lead variant      | Chr:pos (GRCh38)        | Pathway              | Why this locus |
|--------|-------------------|--------------------------|----------------------|----------------|
| CFH    | rs1061170 (Y402H) | 1:196,690,030            | Complement cascade   | Strongest AMD GWAS hit; missense in complement factor H. Cell-intrinsic and circulating complement regulator. |
| ARMS2/HTRA1 | rs10490924 (A69S) | 10:122,454,932      | Extracellular matrix / serine protease | Second-strongest AMD GWAS hit; in LD with the HTRA1 promoter variant rs11200638. |
| C3     | rs2230199 (R102G) | 19:6,718,387             | Complement cascade   | Coding variant in central complement component; AMD-associated. |

## Anchor cell types

- **Retinal pigment epithelium (RPE)** — CL:0002586. The pigmented monolayer behind the
  photoreceptors. Site of drusen accumulation in early/dry AMD. The functional
  consequence of complement dysregulation in AMD is generally agreed to be on RPE.
- **Choroidal endothelial cells** — CL:0002145. The vascular bed underlying RPE/Bruch's
  membrane; relevant to neovascular ("wet") AMD.

## Anchor scRNA atlases

- **Voigt et al. 2019** (Exp Eye Res 184:234) — human RPE/choroid scRNA-seq, healthy.
  Used to define baseline cell-type-resolved gene programs.
- **Orozco et al. 2023** (Cell Rep 42:112243) — single-cell atlas of human retina with
  AMD and control donors. Source of differential expression evidence.
- **Menon et al. 2019** (Nat Commun 10:4902) — additional healthy retina atlas for
  cross-validation.

## The hypothesis the agent should produce

After running the full workflow on objective "analyze AMD GWAS and generate cellular
hypotheses", the agent should output something like:

> **Hypothesis.** The AMD risk variant **rs1061170** (Y402H, PIP ≈ 0.9) drives disease
> via impaired regulation of the **alternative complement pathway** in **retinal
> pigment epithelial cells**. The variant colocalizes with a cis-eQTL/pQTL for **CFH**;
> CFH is a high-confidence member of the Complement Cascade (Reactome R-HSA-166658).
> A scRNA gene program co-expressing complement regulators (CFH, CFI, CFB, C3) is
> active in RPE; that program is upregulated in AMD donors vs controls.
>
> **Evidence chain:**
> 1. rs1061170 (PIP=0.91, GCST003473 / OTG)
> 2. cis-eQTL coloc with CFH in retina (PP4=0.93, eQTL Catalogue r6)
> 3. CFH ∈ Complement Cascade (Reactome R-HSA-166658, evidence=0.98)
> 4. RPE complement gene program (CFH, CFI, CFB, C3; coherence=0.86, Voigt 2019)
> 5. Program upregulated in AMD RPE (mean log2FC=1.4, padj<0.001, Orozco 2023)

If the agent produces this — or an equivalent multi-step evidence chain anchored on
ARMS2/HTRA1 — the architecture is validated.

## Out-of-scope for v0

- WNT pathway story (illustrative in the design doc; not the dominant AMD biology)
- Newer/larger AMD GWAS (Han et al. 2020, etc.)
- Wet-AMD-specific neovascularization story
- Cross-tissue coloc beyond retina/RPE/choroid
- Multiple-objective workflows (only AMD hypothesis generation in v0)
