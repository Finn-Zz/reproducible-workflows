#!/usr/bin/env python3
"""Score the packaged 432 CIFs with the checked-in QMOF model.

Run this script from any directory. It uses the relative paths in the
per-metal manifests and writes a compact report suitable for comparison with
the checked-in score tables.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parents[1]
MODEL = HERE / "model" / "qmof_geometry_v1.json.gz"
METALS = ("TiCo", "TiMg", "TiNi")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metal", choices=(*METALS, "all"), default="all")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    sys.path.insert(0, str(HERE / "src"))
    from mofrelax import ReferenceModel, score_cif  # noqa: PLC0415

    model = ReferenceModel.load(MODEL)
    metals = METALS if args.metal == "all" else (args.metal,)
    reports: list[dict[str, object]] = []
    for metal in metals:
        manifest = HERE / "manifests" / f"{metal}_candidate_manifest.csv"
        for row in read_rows(manifest):
            cif = PACKAGE / row["package_structure_path"]
            report = score_cif(cif, model)
            reports.append(
                {
                    "current_basename": row["current_basename"],
                    "raspa_basename": row["raspa_basename"],
                    "metal": metal,
                    "status": report.get("status", ""),
                    "score": report.get("score", ""),
                    "threshold": report.get("threshold", model.threshold),
                    "coverage": report.get("coverage", ""),
                    "package_structure_path": row["package_structure_path"],
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "current_basename", "raspa_basename", "metal", "status", "score",
        "threshold", "coverage", "package_structure_path",
    ]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(reports)
    print(f"Wrote {len(reports)} reports to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
