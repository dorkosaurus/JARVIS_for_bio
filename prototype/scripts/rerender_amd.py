"""Re-render the 6 AMD targets with the new two-view (domain + protein) PNGs.

Reads prototype/cache/amd_targets.json, calls render_variant_png for each
target with residue + variant label, and writes back the new render dict.
ESM3 fold is cached on disk so this only re-runs PyMOL.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path("/home/ubuntu/JARVIS_for_bio")
sys.path.insert(0, str(REPO))

from prototype.mcp_servers.esm3_server import render_variant_png

EVIDENCE_JSON = REPO / "prototype" / "cache" / "amd_targets.json"


def main() -> int:
    data = json.loads(EVIDENCE_JSON.read_text())
    for sym, ev in data.items():
        vc = ev.get("consequence") or {}
        uniprot = ev["uniprot"]
        if vc.get("found") and vc.get("protein_start"):
            resi = vc["protein_start"]
            aa = vc.get("amino_acids") or ""
            label = None
            if "/" in aa:
                ref, alt = aa.split("/")
                label = f"{ref}{resi}{alt}"
            print(f"[{sym}] {uniprot} r{resi} {label} (two views)")
            result = render_variant_png(
                uniprot_id=uniprot,
                residue_number=resi,
                variant_label=label,
            )
        else:
            print(f"[{sym}] {uniprot} (full protein, non-coding)")
            result = render_variant_png(uniprot_id=uniprot, residue_number=None)
        if "error" in result:
            print(f"  ERROR: {result['error']}")
            continue
        ev["render"] = result
        for k in ("png_path_domain", "png_path_protein", "png_path"):
            if result.get(k):
                print(f"  {k} = {result[k]}")
    EVIDENCE_JSON.write_text(json.dumps(data, indent=2))
    print(f"\nwrote {EVIDENCE_JSON}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
