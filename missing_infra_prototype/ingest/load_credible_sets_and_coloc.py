"""Ingest credible sets + QTL colocalisation from Open Targets Platform GraphQL.

OT Platform (v26.03) merges what used to be OT Genetics. Every AMD GWAS study
that has been fine-mapped (PICS or SuSiE) is exposed as a list of credible
sets, each with:
  - a lead variant + the full locus.rows (variant, PIP, is95CredibleSet)
  - L2G predictions (gene assignments with confidence scores)
  - colocalisation rows against molecular QTLs (h4 = P(shared causal variant))

For v0 we pull the top AMD studies by credible-set count, then insert:
  - all 95%-CS variants into credible_sets, locus label = top L2G symbol
  - all QTL coloc rows (eQTL/pQTL/sQTL/sceQTL) into coloc_results

API: https://api.platform.opentargets.org/api/v4/graphql
"""

from __future__ import annotations

import sqlite3
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

from .init_db import DB_PATH

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

# Pin the OT data release for provenance. Bump when re-ingesting.
OT_RELEASE = "26.03"
SOURCE_DATASET_CS = f"OpenTargets_Platform_{OT_RELEASE}_credibleSets"
SOURCE_DATASET_COLOC = f"OpenTargets_Platform_{OT_RELEASE}_colocalisation"

# AMD studies to pull. Prioritized by credible-set count and field reference status.
AMD_STUDIES = [
    # Fritsche 2016 — the canonical published AMD GWAS, 16k cases
    "GCST003219",
    # FinnGen R12 — largest by case count with full fine-mapping
    "FINNGEN_R12_H7_AMD",
    # Verma 2024 — newer large multi-ancestry, 12k cases
    "GCST90480048",
]

QTL_STUDY_TYPES = ["eqtl", "pqtl", "sqtl", "sceqtl", "scpqtl", "scsqtl"]

PAGE_SIZE = 50
COLOC_PAGE_SIZE = 100
COLOC_PP4_MIN = 0.5  # don't bother storing low-PP4 coloc; trims noise.

# GraphQL query: one credible set with locus + L2G + QTL coloc in one round-trip.
CS_QUERY = """
query($studyIds: [String!], $page: Pagination) {
  credibleSets(studyIds: $studyIds, page: $page) {
    count
    rows {
      studyLocusId
      studyId
      chromosome
      position
      finemappingMethod
      confidence
      pValueMantissa
      pValueExponent
      variant {
        id
        rsIds
        chromosome
        position
        referenceAllele
        alternateAllele
      }
      locus {
        count
        rows {
          variant { id rsIds chromosome position referenceAllele alternateAllele }
          posteriorProbability
          is95CredibleSet
        }
      }
      l2GPredictions {
        rows {
          score
          target { id approvedSymbol }
        }
      }
      study {
        traitFromSource
        publicationFirstAuthor
        publicationDate
        pubmedId
        nCases
        nSamples
      }
    }
  }
}
"""

# Coloc fetched in a follow-up call so we can filter to QTL types only.
COLOC_QUERY = """
query($studyLocusId: String!, $page: Pagination) {
  credibleSet(studyLocusId: $studyLocusId) {
    colocalisation(studyTypes: [eqtl, pqtl, sqtl, sceqtl], page: $page) {
      count
      rows {
        h3
        h4
        colocalisationMethod
        numberColocalisingVariants
        otherStudyLocus {
          studyId
          studyType
          qtlGeneId
          variant { id rsIds }
          study {
            studyType
            projectId
            biosample { biosampleId biosampleName }
            target { id approvedSymbol }
          }
        }
      }
    }
  }
}
"""


def gql(query: str, variables: dict) -> dict:
    for attempt in range(3):
        try:
            resp = requests.post(
                OT_GRAPHQL,
                json={"query": query, "variables": variables},
                timeout=60,
            )
            resp.raise_for_status()
            payload = resp.json()
            if "errors" in payload:
                raise RuntimeError(f"GraphQL errors: {payload['errors']}")
            return payload["data"]
        except Exception as e:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


def fetch_credible_sets(study_id: str) -> list[dict]:
    """Paginate through every credible set for a study."""
    rows: list[dict] = []
    index = 0
    while True:
        data = gql(
            CS_QUERY,
            {"studyIds": [study_id], "page": {"index": index, "size": PAGE_SIZE}},
        )
        batch = data["credibleSets"]["rows"]
        rows.extend(batch)
        total = data["credibleSets"]["count"]
        if len(rows) >= total or not batch:
            break
        index += 1
    return rows


def fetch_coloc(study_locus_id: str) -> list[dict]:
    rows: list[dict] = []
    index = 0
    while True:
        data = gql(
            COLOC_QUERY,
            {
                "studyLocusId": study_locus_id,
                "page": {"index": index, "size": COLOC_PAGE_SIZE},
            },
        )
        batch = data["credibleSet"]["colocalisation"]["rows"]
        rows.extend(batch)
        total = data["credibleSet"]["colocalisation"]["count"]
        if len(rows) >= total or not batch:
            break
        index += 1
    return rows


def locus_label_for(cs: dict) -> str:
    """Use top-scoring L2G gene as the locus label; fall back to chr:pos."""
    l2g = cs.get("l2GPredictions", {}).get("rows", []) or []
    if l2g:
        top = max(l2g, key=lambda r: r["score"])
        return top["target"]["approvedSymbol"]
    v = cs["variant"]
    return f"chr{v['chromosome']}:{v['position']}"


def load(db_path: Path = DB_PATH, studies: list[str] = AMD_STUDIES) -> dict[str, int]:
    now = datetime.now(timezone.utc).isoformat()
    cs_rows: list[tuple] = []
    coloc_rows: list[tuple] = []
    studies_done = 0
    coloc_calls = 0

    for study_id in studies:
        sets = fetch_credible_sets(study_id)
        for cs in sets:
            trait = cs["study"]["traitFromSource"]
            locus = locus_label_for(cs)
            method = cs["finemappingMethod"] or "unknown"

            # Insert one row per variant in the 95% credible set.
            # Some credible sets (e.g. legacy PICS) don't populate locus.rows;
            # fall back to the index variant itself with PIP=1.0.
            locus_block = cs.get("locus") or {}
            locus_rows = locus_block.get("rows") or []
            if not locus_rows:
                v = cs["variant"]
                locus_rows = [
                    {
                        "variant": v,
                        "posteriorProbability": 1.0,
                        "is95CredibleSet": True,
                    }
                ]
            for lr in locus_rows:
                if not lr.get("is95CredibleSet"):
                    continue
                v = lr["variant"]
                rsid = (v["rsIds"] or [None])[0] or v["id"]
                cs_rows.append(
                    (
                        rsid,
                        study_id,
                        trait,
                        locus,
                        float(lr["posteriorProbability"] or 0.0),
                        v["chromosome"],
                        int(v["position"]),
                        v["referenceAllele"],
                        v["alternateAllele"],
                        method,
                        SOURCE_DATASET_CS,
                        now,
                    )
                )

            # Pull QTL coloc for this credible set.
            colocs = fetch_coloc(cs["studyLocusId"])
            coloc_calls += 1
            lead_v = cs["variant"]
            lead_rsid = (lead_v["rsIds"] or [None])[0] or lead_v["id"]
            for c in colocs:
                if (c.get("h4") or 0) < COLOC_PP4_MIN:
                    continue
                other = c["otherStudyLocus"]
                other_study = other.get("study") or {}
                target = other_study.get("target")
                if not target:
                    continue
                biosample = other_study.get("biosample") or {}
                tissue_id = biosample.get("biosampleId") or "UNKNOWN"
                tissue_name = biosample.get("biosampleName") or "unknown"
                qtl_type = other.get("studyType") or other_study.get("studyType") or "qtl"
                coloc_rows.append(
                    (
                        lead_rsid,
                        target["approvedSymbol"],
                        target["id"],
                        tissue_name,
                        tissue_id,
                        qtl_type,
                        float(c["h4"]),
                        f'{other_study.get("projectId") or "OT"}/{other["studyId"]}',
                        SOURCE_DATASET_COLOC,
                        now,
                    )
                )

        studies_done += 1

    with sqlite3.connect(db_path) as conn:
        conn.execute("DELETE FROM credible_sets")
        conn.execute("DELETE FROM coloc_results")
        conn.executemany(
            "INSERT OR REPLACE INTO credible_sets "
            "(variant_id, gwas_id, trait, locus, pip, chrom, pos, ref, alt, "
            " fine_mapper, source_dataset, retrieved_at) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            cs_rows,
        )
        conn.executemany(
            "INSERT OR REPLACE INTO coloc_results "
            "(variant_id, gene_symbol, gene_ensembl, tissue_label, tissue_ontology, "
            " qtl_type, pp4, qtl_dataset, source_dataset, retrieved_at) "
            "VALUES (?,?,?,?,?,?,?,?,?,?)",
            coloc_rows,
        )
        conn.commit()

    return {
        "studies": studies_done,
        "coloc_api_calls": coloc_calls,
        "credible_set_rows": len(cs_rows),
        "coloc_rows": len(coloc_rows),
    }


if __name__ == "__main__":
    stats = load()
    for k, v in stats.items():
        print(f"  {k}: {v}")
