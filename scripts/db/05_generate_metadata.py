#!/usr/bin/env python3
"""Generate mycelium-standard metadata for the Copper colony-measurement dataset.

Reads the populated DuckDB (db/rhodotorula_phenotypes.duckdb) and a sample
preprocessed Parquet to emit:

    data/metadata/copper-colony-measurements/schema.yaml
    data/metadata/copper-colony-measurements/summary_stats.md

Stats are computed from the database (which holds ALL imported rows), not from
a file sample, so counts/missingness are exact. Column descriptions are
grouped by feature family (Shape_*, Intensity_*, TextureGray_*, Color*,
Bbox_*, Grid_*, Metadata_*) so a 180-column export stays readable.

Usage:
    pixi run python 05_generate_metadata.py
"""

import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.db import get_connection  # noqa: E402
from lib.imagename import parse_image_name  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
META_DIR = REPO_ROOT / "data" / "metadata" / "copper-colony-measurements"
GLOB_ROOT = REPO_ROOT / "data" / "preprocessed"

# Columns the DB adds on top of the 180-source-column parquet export.
DB_ADDED_COLUMNS = {"image_name": "VARCHAR NOT NULL REFERENCES image(image_name)",
                    "object_label": "INTEGER NOT NULL", "well_position": "INTEGER"}

# Family -> (description, units, example regex anchor).
FAMILIES = [
    ("Metadata_", "Image/run metadata stamped onto every object by the segmentation export.",
     "mixed", "Metadata_FileSuffix"),
    ("ObjectLabel", "Sequential label of the detected object (colony) within its image.",
     "count", "ObjectLabel"),
    ("Bbox_", "Bounding-box geometry of the object in pixel coordinates.",
     "Pixel (row/col index)", "Bbox_CenterRR"),
    ("Grid_", "Row/column grid position the object was laid out in.",
     "1-based grid index", "Grid_RowNum"),
    ("Shape_", "Morphometric properties of the object (MorphoLibJ).",
     "Pixels (lengths here), pixels^2 (areas); circularity/roundness dimensionless",
     "Shape_Area"),
    ("Intensity_", "Pixel-intensity statistics of the object in the (gray) channel.",
     "Arbitrary intensity units (AU); dimensionless for ratios/CV",
     "Intensity_IntegratedIntensity"),
    ("TextureGray_", "Haralick co-occurrence texture metrics computed from the "
     "gray-level co-occurrence matrix at four orientations and their average (MorphoLibJ).",
     "Dimensionless", "TextureGray_Contrast-avg-scale05"),
    ("Colorxy_", "Median-filtered chromaticity (x,y) statistics of the object.",
     "Dimensionless chromaticity coordinates", "Colorxy_xMean"),
    ("ColorLab_", "Perceptual CIE L*a*b* color statistics of the object.",
     "CIELAB: L* 0-100, a*/b* signed", "ColorLab_a*Mean"),
    ("ColorHSV_", "HSV color statistics of the object (Hue, Saturation, Brightness).",
     "Hue 0-1 (or 0-360); Saturation/Brightness 0-1", "ColorHSV_HueMean"),
]

FAMILY_DESC = dict((f[0], f[1]) for f in FAMILIES)
FAMILY_UNIT = dict((f[0], f[2]) for f in FAMILIES)


def family_of(col: str) -> str:
    for prefix, _, _, _ in FAMILIES:
        if col.startswith(prefix):
            return prefix
    return "_other"


def get_source_schema(con) -> list[tuple[str, str]]:
    """Return (name, parquet_dtype) for the source parquet columns (one representative file)."""
    file = sorted(GLOB_ROOT.glob("*.parquet"))[0]
    rows = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{file}')").fetchall()
    return [(r[0], r[1]) for r in rows]


def column_overview(con) -> dict:
    """Aggregate stats from the DB colony_measurement + image tables."""
    info = {}
    n_imgs, n_measure = con.execute(
        "SELECT (SELECT count(*) FROM image), (SELECT count(*) FROM colony_measurement)"
    ).fetchone()
    info["n_images"] = int(n_imgs)
    info["n_measurements"] = int(n_measure)
    info["n_well_positions"] = int(con.execute(
        "SELECT count(DISTINCT well_position) FROM colony_measurement").fetchone()[0])
    info["n_db_columns"] = int(con.execute("SELECT count(*) FROM information_schema.columns "
                                           "WHERE table_name='colony_measurement'").fetchone()[0])
    info["n_source_columns"] = int(sum(1 for _, _ in get_source_schema(con)))
    byte_size = con.execute(
        "SELECT coalesce(sum(file_size_bytes),0) FROM image").fetchone()[0]
    info["size_bytes"] = int(byte_size)
    info["runs"] = [int(r[0]) for r in con.execute(
        "SELECT run_number FROM imager_run ORDER BY run_number").fetchall()]
    info["per_run_images"] = dict(
        con.execute("SELECT run_number, count(*) FROM image GROUP BY run_number ORDER BY run_number").fetchall())
    info["per_run_measurements"] = dict(
        con.execute("SELECT i.run_number, count(*) FROM colony_measurement cm "
                    "JOIN image i USING (image_name) "
                    "GROUP BY i.run_number ORDER BY i.run_number").fetchall())
    lo, hi = con.execute("SELECT min(imaged_at), max(imaged_at) FROM image").fetchone()
    info["date_min"], info["date_max"] = str(lo)[:10], str(hi)[:10]
    return info


def missingness(con, cols: list[tuple[str, str]]) -> list[tuple[str, int, float]]:
    """Null counts for representative numeric columns, computed in the DB."""
    out = []
    numeric = [(c, dt) for c, dt in cols if dt in ("double", "int64", "uint16", "BIGINT")]
    # Cap the probe to a representative subset so the report stays readable.
    for c, _ in numeric[:40]:
        n = int(con.execute(f'SELECT count(*) FROM colony_measurement WHERE "{c}" IS NULL').fetchone()[0])
        if n > 0:
            out.append((c, n, 100.0 * n / AVAIL["n_measurements"]))
    return out


def top_values(con, cols: list[tuple[str, str]], n: int = 3) -> dict[str, list[str]]:
    top = {}
    names = {c for c, _ in cols}
    for c in ("Metadata_FileSuffix", "Metadata_ImageType", "Grid_RowNum", "Grid_ColNum"):
        if c in names:
            vals = [f"{v} ({cnt})" for v, cnt in con.execute(
                f'SELECT "{c}", count(*) FROM colony_measurement GROUP BY "{c}" '
                f'ORDER BY count(*) DESC LIMIT {n}').fetchall()]
            top[c] = vals
    return top


def write_schema(con, cols: list[tuple[str, str]], meta: dict) -> None:
    lines = [
        "# Dataset Schema", "# Computed by scripts/db/05_generate_metadata.py",
        f"# Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d')}",
        "",
        "dataset_name: copper-colony-measurements",
        "description: >",
        "  Per-colony measurements for the Copper exposure phenotyping experiment (runs 353-357).",
        "  One row per segmented colony object in each imaged plate well. Each of the 2,398",
        "  preprocessed Parquet files holds the objects detected in one image (one run/plate",
        "  timepoint); rows are keyed by image_name + object_label. Source columns (180) plus",
        "  three database-added key columns (image_name, object_label, well_position).",
        "format: Parquet (source), DuckDB colony_measurement table (database)",
        "unit_of_observation: one row per colony object per image",
        "",
        "column_families:",
    ]
    for prefix, desc, units, _ in FAMILIES:
        lines += [f"  {prefix}:", f"    description: {desc}", f"    units: {units}"]
    lines += ["", "columns:"]
    for name, dtype in cols:
        fam = family_of(name)
        lines += [
            f"  - name: \"{name}\"",
            f"    family: {fam}",
            f"    type: \"{dtype}\"",
            f"    description: >",
            f"      {FAMILY_DESC.get(fam, 'Open feature')} See column_families.{fam}.",
            f"    units: {FAMILY_UNIT.get(fam, 'mixed')}",
            f"    nullable: true",
        ]
    lines += ["", "# Database-added key columns (not in the source parquet):"]
    for name, dtype in DB_ADDED_COLUMNS.items():
        lines += [
            f"  - name: \"{name}\"",
            f"    family: _database_key",
            f"    type: \"{dtype}\"",
            f"    description: DB-added link to the image row / object identity and well position.",
            f"    units: n/a",
            f"    nullable: {'true' if name == 'well_position' else 'false'}",
        ]
    META_DIR.mkdir(parents=True, exist_ok=True)
    (META_DIR / "schema.yaml").write_text("\n".join(lines) + "\n")


def write_summary(con, cols: list[tuple[str, str]], meta: dict, miss) -> None:
    def fmt(n: int) -> str:
        return f"{n:,}"
    per_run_rows = "\n".join(
        f"| {r} | {fmt(meta['per_run_images'].get(r,0))} | {fmt(meta['per_run_measurements'].get(r,0))} |"
        for r in meta["runs"])
    freq = top_values(con, cols)
    freq_lines = "\n".join(
        f"- **{c}**: {', '.join(vals)}" for c, vals in freq.items())
    if miss:
        miss_rows = "\n".join(
            f"| {c} | {n:,} | {pct:.2f}% | Computed over all imported rows |" for c, n, pct in miss)
    else:
        miss_rows = "| None | 0 | 0.00% | No numeric probe columns with nulls |"
    n_dropped = 2  # Metadata_Dataset, Metadata_ImageName (see 20_import_measurements.py)
    lines = [
        "# Summary Statistics: copper-colony-measurements",
        "",
        f"Computed by scripts/db/05_generate_metadata.py from DuckDB colony_measurement (all imported rows).",
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%d')}",
        "",
        "## Overview",
        "",
        "| Property | Value |",
        "|----------|-------|",
        f"| Image files (Parquet) | {fmt(meta['n_images'])} |",
        f"| Measurement rows | {fmt(meta['n_measurements'])} |",
        f"| Database columns | {meta['n_db_columns']} ({meta['n_source_columns']} parquet - {n_dropped} dropped + 3 added) |",
        f"| Distinct well positions | {fmt(meta['n_well_positions'])} |",
        f"| Imager runs | {len(meta['runs'])} (d000353-d000357) |",
        f"| Date range | {meta['date_min']} to {meta['date_max']} |",
        f"| Aggregate source size | {fmt(meta['size_bytes'])} bytes |",
        "| Format | Parquet -> DuckDB `colony_measurement` |",
        "",
        "## Per-run coverage",
        "",
        "| Run | Images | Measurement rows |",
        "|-----|--------|------------------|",
        per_run_rows,
        "",
        "## Column family layout",
        "",
        "| Family | Meaning |",
        "|--------|---------|",
    ]
    for prefix, desc, _, _ in FAMILIES:
        cnt = sum(1 for c, _ in cols if family_of(c) == prefix)
        lines.append(f"| `{prefix}` | {cnt} col(s) — {FAMILY_DESC[prefix]} |")
    lines += [
        "",
        "## Missing data summary",
        "",
        "Probe of representative numeric columns (computed over all imported rows):",
        "",
        "| Column | Missing count | Missing % | Pattern / notes |",
        "|--------|---------------|-----------|------------------|",
        miss_rows,
        "",
        "## Most frequent categorical values",
        "",
        freq_lines or "- None probed.",
        "",
        "## Quality flags",
        "",
        "- **Object count can vary by image/timepoint** — wells not always fully detected;",
        "  per-image object counts differ (see variable file sizes in preprocessed/).",
        "- Metadata_Dataset / Metadata_ImageName source columns are dropped on import",
        "  (per-dataset temp identifiers / already encoded in image_name). See DATABASE_DESIGN.md.",
        "",
        "## Notes",
        "",
        "- Image filenames encode run_number, temperature_token, plate_number, imaged_at"
        "  (scripts/db/lib/imagename.py).",
        "- well_position is derived at import from Grid row/col (scripts/db/20_import_measurements.py).",
        "",
    ]
    (META_DIR / "summary_stats.md").write_text("\n".join(lines) + "\n")


def main():
    con = get_connection()
    cols = get_source_schema(con)
    global AVAIL
    AVAIL = column_overview(con)
    miss = missingness(con, cols)
    write_schema(con, cols, AVAIL)
    write_summary(con, cols, AVAIL, miss)
    sizes = list(GLOB_ROOT.glob("*.parquet"))
    print(f"Wrote metadata for {len(sizes)} parquet files under {META_DIR}")
    print(f"measurements={AVAIL['n_measurements']:,} rows across {AVAIL['n_images']:,} images; "
          f"runs={AVAIL['runs']}; date range {AVAIL['date_min']}..{AVAIL['date_max']}")


AVAIL = {}
if __name__ == "__main__":
    main()