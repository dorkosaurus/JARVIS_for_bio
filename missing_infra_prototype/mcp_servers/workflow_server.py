"""JARVIS workflow MCP server.

Exposes the analytic-flow DAG so the agent doesn't have to invent the workflow.
Given an objective the agent gets back an ordered sequence of steps with
explicit tool names, input/output shapes, and dependencies.

Run via:
    ~/venv/bin/python -m prototype.mcp_servers.workflow_server
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, Optional

from mcp.server.fastmcp import FastMCP

DB_PATH = (
    Path(__file__).resolve().parent.parent / "indices" / "jarvis.sqlite"
)

mcp = FastMCP("jarvis-workflow")


def _conn() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


@mcp.tool()
def list_workflows() -> dict[str, Any]:
    """List all analytic workflows available in the workflow DB."""
    with _conn() as c:
        rows = [dict(r) for r in c.execute("SELECT workflow_id, objective, description FROM workflows")]
    return {"count": len(rows), "workflows": rows}


@mcp.tool()
def get_workflow(workflow_id: str) -> dict[str, Any]:
    """Get the full DAG of an analytic workflow.

    Args:
        workflow_id: e.g. "Run_GWAS_Discovery".
    Returns:
        {workflow_id, objective, description, steps: [{step_id, step_order,
         step_name, description, tool_name, input_schema, output_schema,
         depends_on}]}
    """
    with _conn() as c:
        wf_row = c.execute(
            "SELECT workflow_id, objective, description FROM workflows WHERE workflow_id = ?",
            [workflow_id],
        ).fetchone()
        if wf_row is None:
            return {"error": f"workflow {workflow_id!r} not found"}
        steps = [
            dict(r)
            for r in c.execute(
                "SELECT step_id, step_order, step_name, description, tool_name, "
                "       input_schema, output_schema, depends_on "
                "FROM workflow_steps WHERE workflow_id = ? ORDER BY step_order",
                [workflow_id],
            )
        ]
    for s in steps:
        s["input_schema"] = json.loads(s["input_schema"])
        s["output_schema"] = json.loads(s["output_schema"])
    return {**dict(wf_row), "steps": steps}


@mcp.tool()
def find_workflow_for_objective(objective_substring: str) -> dict[str, Any]:
    """Find a workflow whose objective text matches a substring.

    The agent uses this to pick a workflow from a free-text user request.
    """
    with _conn() as c:
        rows = [
            dict(r)
            for r in c.execute(
                "SELECT workflow_id, objective, description FROM workflows "
                "WHERE LOWER(objective) LIKE ? OR LOWER(description) LIKE ?",
                [f"%{objective_substring.lower()}%", f"%{objective_substring.lower()}%"],
            )
        ]
    return {"count": len(rows), "workflows": rows}


if __name__ == "__main__":
    mcp.run()
