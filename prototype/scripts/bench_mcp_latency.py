"""Clean latency benchmark: where does time actually go in v0?

Methodology:
  - 1 warm-up call (untimed) per operation
  - 50 iterations of each operation in steady state
  - Report median + p95 (latency tail matters more than mean)
  - Separate framework overhead (FastMCP dispatch) from underlying work
"""
from __future__ import annotations

import statistics
import subprocess
import sys
import time
from pathlib import Path

import requests

sys.path.insert(0, "/home/ubuntu/JARVIS_for_bio")
from prototype.mcp_servers import ot_server as ot


def percentile(xs: list[float], p: float) -> float:
    xs = sorted(xs)
    k = (len(xs) - 1) * p
    f = int(k)
    c = min(f + 1, len(xs) - 1)
    if f == c:
        return xs[f]
    return xs[f] + (xs[c] - xs[f]) * (k - f)


def time_n(fn, n: int = 50) -> tuple[float, float, float]:
    fn()  # warm-up
    samples = []
    for _ in range(n):
        t = time.perf_counter()
        fn()
        samples.append((time.perf_counter() - t) * 1000)
    return statistics.median(samples), percentile(samples, 0.95), max(samples)


def main() -> int:
    # warm DuckDB attach
    _ = ot._db()

    print(f"{'operation':55s}  {'median':>9s}  {'p95':>9s}  {'max':>9s}")
    print("-" * 90)

    # ---- jarvis-ot, indexed slim cache (the hot path) -------------------
    m, p95, mx = time_n(lambda: ot.study_lookup("GCST003219"))
    print(f"{'ot.study_lookup (indexed cache)':55s}  {m:7.2f} ms  {p95:7.2f} ms  {mx:7.2f} ms")

    m, p95, mx = time_n(lambda: ot.l2g_top_genes("GCST003219", 0.5, 20))
    print(f"{'ot.l2g_top_genes (indexed slim join)':55s}  {m:7.2f} ms  {p95:7.2f} ms  {mx:7.2f} ms")

    m, p95, mx = time_n(lambda: ot.credible_sets_for_study("GCST003219"))
    print(f"{'ot.credible_sets_for_study (indexed cache)':55s}  {m:7.2f} ms  {p95:7.2f} ms  {mx:7.2f} ms")

    # gene_metadata reads target_full (parquet view) — slower
    m, p95, mx = time_n(lambda: ot.gene_metadata("C3"))
    print(f"{'ot.gene_metadata (parquet — transcripts struct)':55s}  {m:7.2f} ms  {p95:7.2f} ms  {mx:7.2f} ms")

    # lead_variant_for_locus reads credible_set_locus (parquet view)
    m, p95, mx = time_n(
        lambda: ot.lead_variant_for_locus("a2fc4eb7d11a5fabbe0e9141a92bcc9a")
    )
    print(f"{'ot.lead_variant_for_locus (parquet — locus struct)':55s}  {m:7.2f} ms  {p95:7.2f} ms  {mx:7.2f} ms")

    # l2g_feature_contributions reads l2g_features (parquet view)
    m, p95, mx = time_n(
        lambda: ot.l2g_feature_contributions(
            "a2fc4eb7d11a5fabbe0e9141a92bcc9a", "ENSG00000125730"
        )
    )
    print(f"{'ot.l2g_feature_contributions (parquet — features)':55s}  {m:7.2f} ms  {p95:7.2f} ms  {mx:7.2f} ms")

    print()

    # ---- FastMCP framework overhead --------------------------------------
    fn = ot.mcp._tool_manager._tools["study_lookup"].fn
    m, p95, mx = time_n(lambda: fn(study_id="GCST003219"))
    print(f"{'FastMCP wrapper around study_lookup':55s}  {m:7.2f} ms  {p95:7.2f} ms  {mx:7.2f} ms")
    print()

    # ---- remote work (sample size = 10, slower) --------------------------
    def t10(fn, n=10):
        fn()
        samples = []
        for _ in range(n):
            t = time.perf_counter()
            fn()
            samples.append((time.perf_counter() - t) * 1000)
        return statistics.median(samples), percentile(samples, 0.95), max(samples)

    print("  (slower remote calls, n=10)")
    print()

    def safe(fn):
        try:
            fn()
            return True
        except Exception as e:
            return False

    def vep_call():
        requests.get(
            "https://rest.ensembl.org/vep/human/id/rs2230199",
            headers={"Accept": "application/json"},
            timeout=60,
        )
    if safe(vep_call):
        m, p95, mx = t10(vep_call)
        print(f"{'Ensembl VEP REST (rs2230199, uncached)':55s}  {m:7.0f} ms  {p95:7.0f} ms  {mx:7.0f} ms")
    else:
        print(f"{'Ensembl VEP REST':55s}  (timeout/rate-limited)")

    def uniprot_call():
        requests.get("https://rest.uniprot.org/uniprotkb/P01024.fasta", timeout=60)
    if safe(uniprot_call):
        m, p95, mx = t10(uniprot_call)
        print(f"{'UniProt REST (P01024 FASTA)':55s}  {m:7.0f} ms  {p95:7.0f} ms  {mx:7.0f} ms")
    else:
        print(f"{'UniProt REST':55s}  (timeout)")

    def paperclip_call():
        subprocess.run(
            [
                "paperclip", "--no-repo", "search",
                "C3 complement AMD", "-n", "3", "-s", "papers",
            ],
            capture_output=True, text=True, timeout=60,
        )
    m, p95, mx = t10(paperclip_call, n=5)
    print(f"{'paperclip search (n=3, source=papers)':55s}  {m:7.0f} ms  {p95:7.0f} ms  {mx:7.0f} ms")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
