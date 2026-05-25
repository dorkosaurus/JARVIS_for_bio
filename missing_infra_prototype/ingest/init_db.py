"""Initialize the SQLite indices database by applying schema.sql."""

import sqlite3
from pathlib import Path

PROTOTYPE_ROOT = Path(__file__).resolve().parent.parent
DB_PATH = PROTOTYPE_ROOT / "indices" / "jarvis.sqlite"
SCHEMA_PATH = PROTOTYPE_ROOT / "ingest" / "schema.sql"


def init_db(db_path: Path = DB_PATH, schema_path: Path = SCHEMA_PATH) -> Path:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    schema_sql = schema_path.read_text()
    with sqlite3.connect(db_path) as conn:
        conn.executescript(schema_sql)
    return db_path


if __name__ == "__main__":
    path = init_db()
    with sqlite3.connect(path) as conn:
        tables = [
            r[0]
            for r in conn.execute(
                "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
            )
        ]
    print(f"Initialized {path}")
    for t in tables:
        print(f"  - {t}")
