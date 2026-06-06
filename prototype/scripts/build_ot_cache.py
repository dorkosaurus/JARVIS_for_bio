"""Build a persistent DuckDB cache over the OT parquets — memory-frugal v2.

Materializes only the slim columns the agent-facing hot path needs:
  - study (slim, no array cols)
  - target (slim, no transcripts/canonicalTranscript struct)
  - study_l2g_slim (the materialized join — (studyId, studyLocusId, position,
    geneId, gene_symbol, l2g_score, biotype))

Struct-heavy columns are intentionally NOT materialized; they stay on parquet:
  - credible_set.locus (variant-membership array) — used by lead_variant_for_locus
  - l2g_prediction.features (29-feature SHAP struct) — used by l2g_feature_contributions
  - target.canonicalTranscript / transcripts — used by gene_metadata

Runs comfortably under 1.5 GB RAM. Output ~150–300 MB at data/ot/cache.duckdb.
"""
from __future__ import annotations

import time
from pathlib import Path

import duckdb

DATA = Path("/home/ubuntu/JARVIS_for_bio/data/ot")
CACHE = DATA / "cache.duckdb"


def main() -> int:
    if CACHE.exists():
        print(f"  removing existing cache: {CACHE}")
        CACHE.unlink()

    con = duckdb.connect(str(CACHE))
    con.execute("PRAGMA memory_limit='1200MB'")
    con.execute("PRAGMA temp_directory='/tmp/duckdb_spill'")
    con.execute("PRAGMA threads=2")
    print(f"  building {CACHE} (memory_limit=1.2GB, spill to /tmp/duckdb_spill)")

    # ---- slim study table ------------------------------------------------
    t = time.perf_counter()
    con.execute(f"""
        CREATE TABLE study AS
        SELECT studyId, traitFromSource, pubmedId,
               publicationFirstAuthor, publicationDate, publicationTitle,
               nCases, nControls, nSamples, studyType, hasSumstats
        FROM read_parquet('{DATA}/study/*.parquet')
    """)
    n = con.execute("SELECT COUNT(*) FROM study").fetchone()[0]
    print(f"    study             : {n:>10,} rows in {time.perf_counter()-t:.1f}s")

    # ---- slim target table (no nested struct arrays) --------------------
    t = time.perf_counter()
    con.execute(f"""
        CREATE TABLE target AS
        SELECT id, approvedSymbol, approvedName, biotype
        FROM read_parquet('{DATA}/target/*.parquet')
    """)
    n = con.execute("SELECT COUNT(*) FROM target").fetchone()[0]
    print(f"    target (slim)     : {n:>10,} rows in {time.perf_counter()-t:.1f}s")

    # ---- materialized study_l2g_slim (the hot path) ---------------------
    # Single-statement insert from streamed parquet read; DuckDB will spill
    # the join hash table if needed.
    print("  materializing study_l2g_slim (streaming join)...")
    t = time.perf_counter()
    con.execute(f"""
        CREATE TABLE study_l2g_slim AS
        SELECT cs.studyId,
               cs.studyLocusId,
               cs.chromosome,
               cs.position,
               l2g.geneId,
               t.approvedSymbol AS gene_symbol,
               t.biotype,
               l2g.score AS l2g_score
        FROM read_parquet('{DATA}/credible_set/*.parquet') cs
        JOIN read_parquet('{DATA}/l2g_prediction/*.parquet') l2g
          ON cs.studyLocusId = l2g.studyLocusId
        LEFT JOIN target t ON l2g.geneId = t.id
    """)
    n = con.execute("SELECT COUNT(*) FROM study_l2g_slim").fetchone()[0]
    print(f"    study_l2g_slim    : {n:>10,} rows in {time.perf_counter()-t:.1f}s")

    # ---- indexes ---------------------------------------------------------
    print("  building indexes...")
    t = time.perf_counter()
    for sql in [
        "CREATE INDEX idx_study_id ON study(studyId)",
        "CREATE INDEX idx_target_id ON target(id)",
        "CREATE INDEX idx_target_symbol ON target(approvedSymbol)",
        "CREATE INDEX idx_sl_studyid ON study_l2g_slim(studyId)",
        "CREATE INDEX idx_sl_slid ON study_l2g_slim(studyLocusId)",
        "CREATE INDEX idx_sl_geneid ON study_l2g_slim(geneId)",
        "CREATE INDEX idx_sl_symbol ON study_l2g_slim(gene_symbol)",
    ]:
        con.execute(sql)
    print(f"    indexes done in {time.perf_counter()-t:.1f}s")

    t = time.perf_counter()
    con.execute("CHECKPOINT")
    print(f"  checkpoint in {time.perf_counter()-t:.1f}s")

    con.close()
    sz = CACHE.stat().st_size / (1024 ** 2)
    print(f"\n  {CACHE}  ({sz:.0f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
