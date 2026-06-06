"""JARVIS PaperClip MCP server (thin wrapper over `paperclip` CLI).

Wraps the GXL paperclip CLI (https://paperclip.gxl.ai) to expose a single
literature-mining tool to the JARVIS workflow. Avoids the multi-step
search→map→fetch race in the raw CLI by stopping at search; the workflow
asks targeted questions per gene, not corpus-wide structured maps.

Run via:
    ~/venv/bin/python -m prototype.mcp_servers.paperclip_server
"""

from __future__ import annotations

import re
import subprocess
from typing import Any, Optional

from mcp.server.fastmcp import FastMCP

mcp = FastMCP("jarvis-paperclip")

PAPER_HEADER_RE = re.compile(
    r"^\s*(\d+)\.\s+(.+)$"
)
META_LINE_RE = re.compile(
    r"^\s*([A-Za-z0-9_]+)\s*·\s*([^·]+?)\s*·\s*([\d-]+)\s*$"
)
DOI_RE = re.compile(r"https?://\S+")
SESSION_RE = re.compile(r"\[s_[a-f0-9]+\]|saved to (s_[a-f0-9]+)")


def _parse_search_output(text: str) -> list[dict[str, Any]]:
    """Parse the paperclip search output into structured papers.

    The CLI emits blocks per paper:
       1. <title>
          <authors>
          <paper_id> · <source> · <YYYY-MM-DD>
          <url or DOI>
          "<summary>"
    """
    # Strip ANSI escapes
    text = re.sub(r"\x1b\[[0-9;]*m", "", text)
    papers: list[dict[str, Any]] = []
    current: Optional[dict[str, Any]] = None
    summary_open = False
    for line in text.splitlines():
        m = PAPER_HEADER_RE.match(line)
        if m:
            if current:
                papers.append(current)
            current = {
                "rank": int(m.group(1)),
                "title": m.group(2).strip(),
                "authors": None,
                "paper_id": None,
                "source": None,
                "date": None,
                "url": None,
                "summary": None,
            }
            summary_open = False
            continue
        if current is None:
            continue
        stripped = line.strip()
        if not stripped:
            continue
        if summary_open:
            current["summary"] = (current["summary"] or "") + " " + stripped.rstrip('"')
            if stripped.endswith('"'):
                summary_open = False
            continue
        meta = META_LINE_RE.match(line)
        if meta:
            current["paper_id"] = meta.group(1).strip()
            current["source"] = meta.group(2).strip()
            current["date"] = meta.group(3).strip()
            continue
        if stripped.startswith(("http://", "https://")):
            current["url"] = stripped
            continue
        if stripped.startswith('"'):
            current["summary"] = stripped.strip('"').strip()
            if not stripped.endswith('"') or len(stripped) < 3:
                summary_open = True
            continue
        if current["authors"] is None:
            current["authors"] = stripped
            continue
        current["authors"] = (current["authors"] or "") + " " + stripped
    if current:
        papers.append(current)
    return papers


@mcp.tool()
def literature_for_gene(
    gene_symbol: str,
    disease_context: str = "",
    n: int = 10,
    source: str = "papers",
) -> dict[str, Any]:
    """Search the biomedical literature for a gene (optionally narrowed by disease).

    Args:
        gene_symbol: HGNC symbol (e.g. "C3", "CFH").
        disease_context: Optional disease phrase (e.g. "age-related macular degeneration").
        n: Max number of papers to return (default 10).
        source: PaperClip source channel. Aggregate "papers" recommended; can
            also pass "pmc", "biorxiv", "medrxiv", "arxiv", "trials/us", "fda".

    Returns parsed papers with title, authors, paper_id, date, url, summary.
    """
    query = gene_symbol
    if disease_context:
        query = f"{gene_symbol} {disease_context}"

    try:
        result = subprocess.run(
            ["paperclip", "--no-repo", "search", query, "-n", str(n), "-s", source],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except subprocess.TimeoutExpired:
        return {"error": "paperclip search timed out (60s)"}

    if result.returncode != 0:
        return {
            "error": "paperclip search failed",
            "stderr": result.stderr[-400:],
            "stdout": result.stdout[-400:],
        }

    papers = _parse_search_output(result.stdout)
    return {
        "gene_symbol": gene_symbol,
        "query": query,
        "source": source,
        "count": len(papers),
        "papers": papers[:n],
        "provenance": "PaperClip (paperclip.gxl.ai) — BM25 + vector search over public scientific corpus",
    }


@mcp.tool()
def literature_for_variant(
    variant_label: str,
    gene_symbol: str = "",
    n: int = 5,
    source: str = "papers",
) -> dict[str, Any]:
    """Search for papers about a specific variant (e.g. "C3 R102C", "rs2230199").

    Args:
        variant_label: Variant identifier or HGVS-style label.
        gene_symbol: Optional gene to disambiguate.
        n: Max papers.
        source: PaperClip source channel.
    """
    query = variant_label if not gene_symbol else f"{gene_symbol} {variant_label}"
    try:
        result = subprocess.run(
            ["paperclip", "--no-repo", "search", query, "-n", str(n), "-s", source],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except subprocess.TimeoutExpired:
        return {"error": "paperclip search timed out (60s)"}
    if result.returncode != 0:
        return {"error": "paperclip search failed", "stderr": result.stderr[-400:]}

    papers = _parse_search_output(result.stdout)
    return {
        "variant_label": variant_label,
        "gene_symbol": gene_symbol,
        "query": query,
        "count": len(papers),
        "papers": papers[:n],
        "provenance": "PaperClip (paperclip.gxl.ai)",
    }


if __name__ == "__main__":
    mcp.run()
