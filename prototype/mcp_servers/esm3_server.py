"""JARVIS ESM3 MCP server.

Tools:
  variant_consequence       — VEP REST lookup for a variant (missense/intronic/...)
  uniprot_sequence          — fetch canonical FASTA
  fold_and_annotate         — ESM3 Forge: structure + function annotations
  render_variant_png        — PyMOL headless render with variant residue highlighted
  score_target              — one-stop orchestration: lookup → fold → render → druggability

Forge calls are slow (~30-60s per protein). Results are cached at
  prototype/cache/esm3/<uniprot>/{structure.pdb, function.json, render_<resi>.png}

Run via:
    ~/venv/bin/python -m prototype.mcp_servers.esm3_server
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any, Optional

import numpy as np
import requests
from dotenv import load_dotenv
from mcp.server.fastmcp import FastMCP

REPO = Path("/home/ubuntu/JARVIS_for_bio")
CACHE = REPO / "prototype" / "cache" / "esm3"
CACHE.mkdir(parents=True, exist_ok=True)

# Forge env: zero_to_one stores it there, v0_release may override
for env_path in [
    REPO / "v0_release" / ".env",
    Path("/home/ubuntu/zero_to_one/esm3/.env"),
]:
    if env_path.exists():
        load_dotenv(env_path)

ESM3_API_KEY = os.environ.get("ESM3_API_KEY") or os.environ.get("ESM_API_KEY")
MODEL = "esm3-open-2024-03"

mcp = FastMCP("jarvis-esm3")


# ---------------------------------------------------------------------------
# Variant consequence (Ensembl VEP REST)
# ---------------------------------------------------------------------------

@mcp.tool()
def variant_consequence(variant_id: str) -> dict[str, Any]:
    """Look up VEP-style consequence for a variant.

    Args:
        variant_id: One of:
          - rsID ("rs2230199")
          - OT-style ("19_6718376_G_C")
          - chr:pos:ref/alt ("19:6718376:G/C")

    Returns most_severe_consequence, gene/transcript/protein context, and for
    missense variants: protein position + amino acid change (e.g. "R/G").
    """
    cache_key = CACHE / "vep" / f"{variant_id.replace('/', '_')}.json"
    cache_key.parent.mkdir(parents=True, exist_ok=True)
    if cache_key.exists():
        return json.loads(cache_key.read_text())

    if variant_id.startswith("rs"):
        url = f"https://rest.ensembl.org/vep/human/id/{variant_id}?canonical=1"
    else:
        # Normalise OT-style "19_6718376_G_C" -> "19 6718376 6718376 G/C"
        if "_" in variant_id:
            parts = variant_id.split("_")
            chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
            region = f"{chrom}:{pos}-{pos}:1/{alt}"
        elif ":" in variant_id and "/" in variant_id:
            chrom, pos, alleles = variant_id.split(":")
            ref, alt = alleles.split("/")
            region = f"{chrom}:{pos}-{pos}:1/{alt}"
        else:
            return {"error": f"unrecognized variant_id format: {variant_id}"}
        url = f"https://rest.ensembl.org/vep/human/region/{region}?canonical=1"

    r = requests.get(url, headers={"Accept": "application/json"}, timeout=30)
    if not r.ok:
        return {"error": f"VEP returned {r.status_code}: {r.text[:200]}"}
    data = r.json()
    if not data:
        return {"found": False, "variant_id": variant_id}
    rec = data[0]

    # Find the canonical transcript consequence
    canon = None
    for tc in rec.get("transcript_consequences", []):
        if tc.get("canonical") == 1:
            canon = tc
            break
    canon = canon or (rec.get("transcript_consequences") or [{}])[0]

    result = {
        "found": True,
        "variant_id": variant_id,
        "most_severe_consequence": rec.get("most_severe_consequence"),
        "assembly": rec.get("assembly_name"),
        "chromosome": str(rec.get("seq_region_name", "")),
        "position": rec.get("start"),
        "ref_allele": rec.get("allele_string", "").split("/")[0] if rec.get("allele_string") else None,
        "alt_allele": rec.get("allele_string", "").split("/")[-1] if rec.get("allele_string") else None,
        "gene_symbol": canon.get("gene_symbol"),
        "gene_id": canon.get("gene_id"),
        "transcript_id": canon.get("transcript_id"),
        "consequence_terms": canon.get("consequence_terms", []),
        "protein_start": canon.get("protein_start"),
        "protein_end": canon.get("protein_end"),
        "amino_acids": canon.get("amino_acids"),
        "codons": canon.get("codons"),
        "biotype": canon.get("biotype"),
        "impact": canon.get("impact"),
        "polyphen_prediction": canon.get("polyphen_prediction"),
        "polyphen_score": canon.get("polyphen_score"),
        "sift_prediction": canon.get("sift_prediction"),
        "sift_score": canon.get("sift_score"),
        "provenance": "Ensembl VEP REST (GRCh38)",
    }
    cache_key.write_text(json.dumps(result))
    return result


# ---------------------------------------------------------------------------
# UniProt sequence
# ---------------------------------------------------------------------------

@mcp.tool()
def uniprot_sequence(uniprot_id: str) -> dict[str, Any]:
    """Fetch canonical protein sequence (FASTA) from UniProt."""
    cache_key = CACHE / "uniprot" / f"{uniprot_id}.json"
    cache_key.parent.mkdir(parents=True, exist_ok=True)
    if cache_key.exists():
        return json.loads(cache_key.read_text())

    r = requests.get(
        f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta",
        timeout=30,
    )
    if not r.ok:
        return {"error": f"UniProt returned {r.status_code}"}
    lines = r.text.strip().split("\n")
    header = lines[0]
    seq = "".join(lines[1:])
    result = {
        "uniprot_id": uniprot_id,
        "header": header,
        "sequence": seq,
        "length": len(seq),
        "provenance": "UniProt REST",
    }
    cache_key.write_text(json.dumps(result))
    return result


# ---------------------------------------------------------------------------
# ESM3 fold + function annotation
# ---------------------------------------------------------------------------

def _esm3_run(sequence: str) -> tuple[str, dict]:
    """Run ESM3 Forge fold + function annotation. Returns (pdb_string, meta)."""
    from esm.sdk import client as make_client
    from esm.sdk.api import ESMProtein, GenerationConfig

    if not ESM3_API_KEY:
        raise RuntimeError("ESM3_API_KEY not set (looked in .env)")
    cli = make_client(model=MODEL, token=ESM3_API_KEY)
    folded = cli.generate(
        ESMProtein(sequence=sequence),
        GenerationConfig(track="structure"),
    )
    if hasattr(folded, "error_code"):
        raise RuntimeError(f"fold error: {folded.error_msg}")
    annotated = cli.generate(
        ESMProtein(sequence=sequence),
        GenerationConfig(track="function"),
    )
    if hasattr(annotated, "error_code"):
        raise RuntimeError(f"function error: {annotated.error_msg}")
    pdb = folded.to_pdb_string()
    plddt = folded.plddt
    if hasattr(plddt, "detach"):
        plddt = plddt.detach().cpu().numpy().tolist()
    elif isinstance(plddt, np.ndarray):
        plddt = plddt.tolist()
    meta = {
        "sequence": sequence,
        "length": len(sequence),
        "plddt": plddt,
        "mean_plddt": float(np.mean(plddt)) if plddt else None,
        "ptm": float(folded.ptm) if folded.ptm is not None else None,
        "function_annotations": [
            {"label": fa.label, "start": fa.start, "end": fa.end}
            for fa in (annotated.function_annotations or [])
        ],
    }
    return pdb, meta


@mcp.tool()
def fold_and_annotate(uniprot_id: str, force: bool = False) -> dict[str, Any]:
    """Run ESM3 fold + InterPro function annotation for a protein.

    Cached. Returns mean pLDDT, pTM, function annotations, and paths to the
    PDB structure and function JSON.
    """
    target_dir = CACHE / uniprot_id
    target_dir.mkdir(exist_ok=True)
    pdb_path = target_dir / "structure.pdb"
    fn_path = target_dir / "function.json"

    if not force and pdb_path.exists() and fn_path.exists():
        meta = json.loads(fn_path.read_text())
        return {
            "uniprot_id": uniprot_id,
            "pdb_path": str(pdb_path),
            "function_path": str(fn_path),
            "length": meta.get("length"),
            "mean_plddt": meta.get("mean_plddt"),
            "ptm": meta.get("ptm"),
            "n_function_annotations": len(meta.get("function_annotations") or []),
            "function_annotations": meta.get("function_annotations"),
            "cached": True,
            "provenance": f"ESM3 Forge ({MODEL})",
        }

    seq_info = uniprot_sequence(uniprot_id)
    if "error" in seq_info:
        return seq_info
    seq = seq_info["sequence"]

    pdb, meta = _esm3_run(seq)
    pdb_path.write_text(pdb)
    fn_path.write_text(json.dumps(meta))
    return {
        "uniprot_id": uniprot_id,
        "pdb_path": str(pdb_path),
        "function_path": str(fn_path),
        "length": meta["length"],
        "mean_plddt": meta["mean_plddt"],
        "ptm": meta["ptm"],
        "n_function_annotations": len(meta["function_annotations"]),
        "function_annotations": meta["function_annotations"],
        "cached": False,
        "provenance": f"ESM3 Forge ({MODEL})",
    }


# ---------------------------------------------------------------------------
# PyMOL headless variant render
# ---------------------------------------------------------------------------

@mcp.tool()
def render_variant_png(
    uniprot_id: str,
    residue_number: Optional[int] = None,
    variant_label: Optional[str] = None,
    width: int = 900,
    height: int = 700,
    zoom_buffer: float = 25.0,
) -> dict[str, Any]:
    """Render the protein structure as cartoon PNG(s).

    If residue_number is given, produces TWO PNGs:
      - render_r<resi>_domain.png  — zoomed on the residue + zoom_buffer Å
        of surrounding context (default 25 Å, i.e. whole local domain)
      - render_r<resi>_protein.png — full protein with the variant residue
        still marked (red sphere + label, sphere scale boosted so it's
        visible at full zoom)
    If residue_number is None, produces one render_all.png of the whole
    protein (no variant marker).

    Outputs: prototype/cache/esm3/<uniprot>/render_*.png
    """
    target_dir = CACHE / uniprot_id
    pdb_path = target_dir / "structure.pdb"
    if not pdb_path.exists():
        return {"error": f"no PDB at {pdb_path} — run fold_and_annotate first"}

    import pymol
    pymol.finish_launching(["pymol", "-cq"])
    from pymol import cmd

    def _setup_scene(sphere_scale: float) -> None:
        cmd.reinitialize()
        cmd.load(str(pdb_path), "prot")
        cmd.bg_color("white")
        cmd.hide("everything")
        cmd.show("cartoon", "prot")
        cmd.color("lightblue", "prot")
        cmd.set("cartoon_transparency", 0.0)
        cmd.set("ray_shadows", 0)
        cmd.set("antialias", 2)
        if residue_number is not None:
            cmd.select("variant", f"resi {residue_number}")
            cmd.show("sticks", "variant")
            cmd.show("spheres", "variant and name CA")
            cmd.set("sphere_scale", sphere_scale, "variant and name CA")
            cmd.color("red", "variant")
            if variant_label:
                cmd.label("variant and name CA", f"'{variant_label}'")
                cmd.set("label_color", "black")
                cmd.set("label_size", 18)

    if residue_number is None:
        png_path = target_dir / "render_all.png"
        _setup_scene(sphere_scale=0.6)
        cmd.zoom("prot")
        cmd.orient("prot")
        cmd.png(str(png_path), width=width, height=height, dpi=150, ray=1)
        return {
            "uniprot_id": uniprot_id,
            "residue_number": None,
            "variant_label": None,
            "png_path": str(png_path),
            "provenance": "PyMOL open-source 3.x headless render",
        }

    # Coding variant: two views
    domain_path = target_dir / f"render_r{residue_number}_domain.png"
    protein_path = target_dir / f"render_r{residue_number}_protein.png"

    # Domain view: orient first to set rotation, then zoom with buffer so
    # the orient call doesn't override the zoom-out.
    _setup_scene(sphere_scale=0.6)
    cmd.orient("variant")
    cmd.zoom("variant", zoom_buffer)
    cmd.png(str(domain_path), width=width, height=height, dpi=150, ray=1)

    # Protein view: big sphere so the variant is visible at full zoom
    _setup_scene(sphere_scale=2.0)
    cmd.orient("prot")
    cmd.zoom("prot")
    cmd.png(str(protein_path), width=width, height=height, dpi=150, ray=1)

    return {
        "uniprot_id": uniprot_id,
        "residue_number": residue_number,
        "variant_label": variant_label,
        "png_path": str(domain_path),  # back-compat
        "png_path_domain": str(domain_path),
        "png_path_protein": str(protein_path),
        "zoom_buffer_angstroms": zoom_buffer,
        "provenance": "PyMOL open-source 3.x headless render",
    }


# ---------------------------------------------------------------------------
# One-stop orchestration
# ---------------------------------------------------------------------------

@mcp.tool()
def score_target(
    gene_symbol: str,
    uniprot_id: str,
    variant_id: Optional[str] = None,
) -> dict[str, Any]:
    """Full evidence pack for one (gene, variant) pair.

    1. VEP lookup of the variant (if given) — finds missense vs non-coding
    2. ESM3 fold + function annotation on the canonical protein
    3. PyMOL render: if variant is coding (has protein_start), zoom on residue
       with the variant labelled; else render the whole protein

    Returns combined dict with variant_consequence, fold summary, and PNG path.
    Use this as the single call per target in the AMD workflow.
    """
    out: dict[str, Any] = {
        "gene_symbol": gene_symbol,
        "uniprot_id": uniprot_id,
        "variant_id": variant_id,
    }

    vc = None
    if variant_id:
        vc = variant_consequence(variant_id)
        out["variant_consequence"] = vc

    fold = fold_and_annotate(uniprot_id)
    out["fold"] = fold
    if "error" in fold:
        return out

    residue = None
    label = None
    if vc and vc.get("found") and vc.get("protein_start"):
        residue = vc["protein_start"]
        aa = vc.get("amino_acids", "")
        if "/" in aa:
            ref_aa, alt_aa = aa.split("/")
            label = f"{ref_aa}{residue}{alt_aa}"

    png = render_variant_png(
        uniprot_id=uniprot_id,
        residue_number=residue,
        variant_label=label,
    )
    out["render"] = png
    out["provenance_chain"] = [
        "Ensembl VEP REST → variant consequence",
        f"UniProt REST → canonical FASTA ({uniprot_id})",
        f"ESM3 Forge ({MODEL}) → fold + function",
        "PyMOL open-source headless → cartoon PNG",
    ]
    return out


if __name__ == "__main__":
    mcp.run()
