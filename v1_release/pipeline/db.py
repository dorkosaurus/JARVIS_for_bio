"""SQLite IO for the closed-loop results.

Schema mirrors CLAUDE.md. The DB is the source of truth; the parquet +
csv exports in `outputs/` are derived from it.
"""

import sqlite3
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config

SCHEMA = """
CREATE TABLE IF NOT EXISTS candidates (
    candidate_id TEXT PRIMARY KEY,
    capsid_id TEXT NOT NULL,
    vp1_sequence TEXT NOT NULL,
    has_7mer_insertion BOOLEAN,
    insertion_peptide TEXT,
    hamming_to_aav2 INTEGER,
    capsid_embedding_path TEXT,
    is_anchor BOOLEAN DEFAULT FALSE,
    anchor_label TEXT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS experimental_results (
    result_id INTEGER PRIMARY KEY AUTOINCREMENT,
    candidate_id TEXT REFERENCES candidates(candidate_id),
    cycle INTEGER NOT NULL,
    selection_strategy TEXT,
    rpe_transduction REAL,
    neut_escape REAL,
    inflammation_score REAL,
    meets_constraint BOOLEAN,
    is_on_pareto_frontier BOOLEAN,
    source TEXT NOT NULL,
    source_version TEXT,
    is_simulated BOOLEAN DEFAULT TRUE,
    assay_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS cycle_summary (
    cycle INTEGER,
    selection_strategy TEXT,
    n_candidates_tested INTEGER,
    n_constraint_violations INTEGER,
    best_rpe_transduction REAL,
    best_neut_escape REAL,
    pareto_frontier_size INTEGER,
    pareto_hypervolume REAL,
    timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    PRIMARY KEY (cycle, selection_strategy)
);
"""


def init_db(path: Path = config.RESULTS_DB, fresh: bool = True) -> sqlite3.Connection:
    """Create or open the SQLite database. If fresh=True, drop existing tables first."""
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(path)
    if fresh:
        conn.executescript(
            "DROP TABLE IF EXISTS experimental_results; "
            "DROP TABLE IF EXISTS cycle_summary; "
            "DROP TABLE IF EXISTS candidates;"
        )
    conn.executescript(SCHEMA)
    conn.commit()
    return conn


def write_candidates(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    """Insert candidate rows. df columns must match the `candidates` table."""
    cols = [
        "candidate_id", "capsid_id", "vp1_sequence",
        "has_7mer_insertion", "insertion_peptide", "hamming_to_aav2",
        "capsid_embedding_path", "is_anchor", "anchor_label",
    ]
    rows = [tuple(r) for r in df[cols].itertuples(index=False, name=None)]
    conn.executemany(
        f"INSERT OR REPLACE INTO candidates ({','.join(cols)}) VALUES ({','.join('?'*len(cols))})",
        rows,
    )
    conn.commit()


def write_results(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    """Insert experimental_results rows."""
    cols = [
        "candidate_id", "cycle", "selection_strategy",
        "rpe_transduction", "neut_escape", "inflammation_score",
        "meets_constraint", "is_on_pareto_frontier",
        "source", "source_version", "is_simulated",
    ]
    rows = [tuple(r) for r in df[cols].itertuples(index=False, name=None)]
    conn.executemany(
        f"INSERT INTO experimental_results ({','.join(cols)}) VALUES ({','.join('?'*len(cols))})",
        rows,
    )
    conn.commit()


def write_cycle_summary(conn: sqlite3.Connection, df: pd.DataFrame) -> None:
    """Insert cycle_summary rows."""
    cols = [
        "cycle", "selection_strategy",
        "n_candidates_tested", "n_constraint_violations",
        "best_rpe_transduction", "best_neut_escape",
        "pareto_frontier_size", "pareto_hypervolume",
    ]
    rows = [tuple(r) for r in df[cols].itertuples(index=False, name=None)]
    conn.executemany(
        f"INSERT OR REPLACE INTO cycle_summary ({','.join(cols)}) VALUES ({','.join('?'*len(cols))})",
        rows,
    )
    conn.commit()


def export_pareto_parquet(conn: sqlite3.Connection, path: Path) -> pd.DataFrame:
    """Join candidates + experimental_results and write the figure-source parquet."""
    df = pd.read_sql_query(
        """
        SELECT
            c.candidate_id, c.capsid_id, c.vp1_sequence,
            c.has_7mer_insertion, c.insertion_peptide, c.hamming_to_aav2,
            c.is_anchor, c.anchor_label,
            r.cycle, r.selection_strategy,
            r.rpe_transduction, r.neut_escape, r.inflammation_score,
            r.meets_constraint, r.is_on_pareto_frontier,
            r.source, r.source_version, r.is_simulated
        FROM candidates c
        JOIN experimental_results r ON c.candidate_id = r.candidate_id
        """,
        conn,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    return df


def export_cycle_summary_csv(conn: sqlite3.Connection, path: Path) -> pd.DataFrame:
    df = pd.read_sql_query("SELECT * FROM cycle_summary ORDER BY selection_strategy, cycle", conn)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return df
