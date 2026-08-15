"""Shared DuckDB connection helper for the import/query scripts."""

from pathlib import Path
import duckdb

DEFAULT_DB_PATH = Path(__file__).resolve().parents[3] / "db" / "rhodotorula_phenotypes.duckdb"


def get_connection(db_path=None) -> duckdb.DuckDBPyConnection:
    db_path = Path(db_path) if db_path else DEFAULT_DB_PATH
    db_path.parent.mkdir(parents=True, exist_ok=True)
    return duckdb.connect(str(db_path))
