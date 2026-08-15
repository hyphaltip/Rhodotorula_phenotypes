"""Parse the Metadata_ImageName filename convention documented in README.md.

Example: d000355_300_079_2026-02-24_17-14-27
  -> run_number=355, temperature_token=300, plate_number=79,
     imaged_at=2026-02-24 17:14:05
"""

import re
from datetime import datetime

IMAGE_NAME_RE = re.compile(
    r"^d(?P<run_number>\d+)_(?P<temperature_token>\d+)_(?P<plate_number>\d+)_"
    r"(?P<datetime>\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2})$"
)


def parse_image_name(image_name: str) -> dict:
    m = IMAGE_NAME_RE.match(image_name)
    if not m:
        raise ValueError(f"Cannot parse Metadata_ImageName: {image_name!r}")
    return {
        "run_number": int(m.group("run_number")),
        "temperature_token": int(m.group("temperature_token")),
        "plate_number": int(m.group("plate_number")),
        "imaged_at": datetime.strptime(m.group("datetime"), "%Y-%m-%d_%H-%M-%S"),
    }
