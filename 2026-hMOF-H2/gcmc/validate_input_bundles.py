#!/usr/bin/env python3
"""Validate the three checked-in RASPA3 example bundles."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
METALS = ("TiCo", "TiMg", "TiNi")


def main() -> int:
    for metal in METALS:
        directory = ROOT / metal
        for filename in ("h2.json", "force_field.json", "simulation.json"):
            path = directory / filename
            if not path.is_file():
                raise FileNotFoundError(path)
        simulation = json.loads((directory / "simulation.json").read_text(encoding="utf-8"))
        system = simulation["Systems"][0]
        cif = directory / f"{system['Name']}.cif"
        if not cif.is_file():
            raise FileNotFoundError(cif)
        if simulation["ForceField"] != "." or simulation["Components"][0]["MoleculeDefinition"] != ".":
            raise ValueError(f"local JSON references changed for {metal}")
        print(f"{metal}: {system['Name']} | {system['ExternalTemperature']} K | {system['ExternalPressure']} Pa")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
