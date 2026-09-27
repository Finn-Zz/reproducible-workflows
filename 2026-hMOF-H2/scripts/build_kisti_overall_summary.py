#!/usr/bin/env python3
"""Build one complete KISTI RASPA summary directly from the ZIP archives.

This bundled module is used as a parser dependency by the v5 rebuild utility;
it is not the standalone package analysis entry point.

Run one temperature at a time to keep memory and runtime bounded, then combine:
  python build_kisti_overall_summary.py --temperature 77
  ...
  python build_kisti_overall_summary.py --combine
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import zipfile
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NEW_ROOT = ROOT.parent / "new"
PARTS = Path("/tmp/kisti_adsorption_summary_parts")
FINAL = ROOT / "KISTI_adsorption_overall_summary.csv"
TEMPERATURES = (77, 120, 160, 200, 233, 253, 273, 298)
ARCHIVE_TEMPERATURES = TEMPERATURES
NUMBER = r"[+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?"
LOADING = re.compile(
    rf"(?P<kind>Abs\.|Excess) loading average\s+(?P<value>{NUMBER})\s+\+/-\s*(?P<error>{NUMBER})\s+\[(?P<unit>molecules/cell|molecules/uc|mol/kg-framework|mg/g-framework)\]"
)
ENTHALPY_K = re.compile(rf"Enthalpy of adsorption:\s*(?P<value>{NUMBER})\s+\+/-\s*(?P<error>{NUMBER})\s+\[K\]")
ENTHALPY_KJ = re.compile(rf"^\s*(?P<value>{NUMBER})\s+\+/-\s*(?P<error>{NUMBER})\s+\[kJ/mol\]", re.MULTILINE)

BASE_FIELDS = (
    "basename", "connectivity", "topology", "metal_node", "organic_linker",
    "temperature_K", "pressure_bar",
    "completed", "completed_via", "result_parse_status",
)
KINDS = ("absolute", "excess")
UNITS = (
    ("molecules_cell", "molecules/cell"),
    ("molecules_uc", "molecules/uc"),
    ("mol_kg_framework", "mol/kg-framework"),
    ("mg_g_framework", "mg/g-framework"),
)
LOADING_FIELDS = tuple(
    f"{kind}_loading_{unit_key}"
    for kind in KINDS
    for unit_key, _ in UNITS
) + (
    "absolute_loading_wt_percent", "excess_loading_wt_percent",
)
EXTRA_FIELDS = (
    "raw_enthalpy_adsorption_K", "raw_enthalpy_adsorption_kJ_mol",
)
FIELDS = BASE_FIELDS + LOADING_FIELDS + EXTRA_FIELDS


def normalize_topology(name: str) -> tuple[str, bool]:
    retry = bool(re.search(r"_(?:incomplete|incompelete)$", name))
    clean = re.sub(r"_(?:incomplete|incompelete)$", "", name).rstrip("_")
    return clean, retry


def expected_jobs() -> list[tuple[str, str, str]]:
    paths = subprocess.check_output(
        ["rg", "--files", "-g", "simulation.json", str(NEW_ROOT)], text=True
    ).splitlines()
    jobs = set()
    for path_text in paths:
        path = Path(path_text)
        topology, _ = normalize_topology(path.parent.parent.parent.name)
        pressure = path.parent.parent.name.rsplit("_", 1)[-1]
        jobs.add((topology, pressure, path.parent.name))
    return sorted(jobs)


def connectivity(basename: str) -> str:
    if "-3c-" in basename:
        return "3,6"
    if "-4c-" in basename:
        return "4,6"
    return "unknown"


def basename_components(basename: str) -> tuple[str, str]:
    match = re.search(r"-(TiCo|TiMg|TiNi)-(N\d+)-", basename)
    return match.groups() if match else ("unknown", "unknown")


def parse_output(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    unit_keys = {unit: key for key, unit in UNITS}
    for match in LOADING.finditer(text):
        kind = "absolute" if match.group("kind") == "Abs." else "excess"
        stem = f"{kind}_loading_{unit_keys[match.group('unit')]}"
        values[stem] = match.group("value")
    for kind in KINDS:
        source = f"{kind}_loading_mg_g_framework"
        if source in values:
            values[f"{kind}_loading_wt_percent"] = f"{float(values[source]) / 10.0:.10g}"
    match = ENTHALPY_K.search(text)
    if match:
        values["raw_enthalpy_adsorption_K"] = match.group("value")
    match = ENTHALPY_KJ.search(text)
    if match:
        values["raw_enthalpy_adsorption_kJ_mol"] = match.group("value")
    return values


def process_temperature(temperature: int, shard_index: int, shard_count: int) -> None:
    if temperature not in ARCHIVE_TEMPERATURES:
        raise SystemExit(f"No ZIP archive is expected for {temperature} K")
    archive = ROOT / f"{temperature}.zip"
    all_expected = expected_jobs()
    expected = [job for index, job in enumerate(all_expected) if index % shard_count == shard_index]
    PARTS.mkdir(parents=True, exist_ok=True)
    success_sources: dict[tuple[str, str, str], list[tuple[bool, str]]] = defaultdict(list)
    output_members: dict[str, str] = {}
    expected_output = {"5": "5e+05", "100": "1e+07"}

    with zipfile.ZipFile(archive) as zf:
        names = zf.namelist()
        for name in names:
            parts = name.rstrip("/").split("/")
            if len(parts) == 4 and parts[-1] == f".raspa_success_{temperature}":
                topology, retry = normalize_topology(parts[0])
                pressure = parts[1].rsplit("_", 1)[-1]
                success_sources[(topology, pressure, parts[2])].append((retry, "/".join(parts[:3])))
            elif len(parts) == 5 and re.fullmatch(rf"output_{temperature}_[0-9]+e\+[0-9]+\.s0\.txt", parts[-1]):
                output_members["/".join(parts[:3])] = name

        rows = []
        for topology, pressure, basename in expected:
            key = (topology, pressure, basename)
            sources = success_sources.get(key, [])
            origin_types = {retry for retry, _ in sources}
            if origin_types == {False, True}:
                completed_via = "first+second"
            elif origin_types == {True}:
                completed_via = "second"
            elif origin_types == {False}:
                completed_via = "first"
            else:
                completed_via = "not_completed"
            row = {field: "" for field in FIELDS}
            metal_node, organic_linker = basename_components(basename)
            row.update({
                "basename": basename, "connectivity": connectivity(basename), "topology": topology,
                "metal_node": metal_node, "organic_linker": organic_linker,
                "temperature_K": str(temperature), "pressure_bar": pressure,
                "completed": "yes" if sources else "no", "completed_via": completed_via,
                "result_parse_status": "not_completed" if not sources else "success_marker_but_output_missing",
            })
            # Prefer a successful second-round output; fall back to first-round.
            chosen = next(((retry, prefix) for retry, prefix in sorted(sources, reverse=True) if prefix in output_members), None)
            if chosen is not None:
                member = output_members[chosen[1]]
                text = zf.read(member).decode("utf-8", errors="replace")
                parsed = parse_output(text)
                row.update(parsed)
                row["result_parse_status"] = "parsed" if "absolute_loading_mg_g_framework" in parsed else "output_unparseable"
            rows.append(row)

    part = PARTS / f"{temperature}.{shard_index}-of-{shard_count}.csv"
    with part.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    completed = sum(row["completed"] == "yes" for row in rows)
    parsed = sum(row["result_parse_status"] == "parsed" for row in rows)
    print(f"{temperature} K shard {shard_index + 1}/{shard_count}: rows={len(rows)}, completed={completed}, parsed={parsed}")


def combine(shard_count: int) -> None:
    expected = expected_jobs()
    rows = []
    for temperature in TEMPERATURES:
        for shard_index in range(shard_count):
            part = PARTS / f"{temperature}.{shard_index}-of-{shard_count}.csv"
            if not part.is_file():
                raise SystemExit(f"Missing temporary result: {part}")
            with part.open(newline="") as handle:
                rows.extend(csv.DictReader(handle))
    # Older temporary parts may predate the descriptive columns and may contain
    # uncertainty columns.  Normalize them here so the public output remains a
    # single compact table with the current schema.
    normalized_rows = []
    for row in rows:
        metal_node, organic_linker = basename_components(row["basename"])
        row["metal_node"] = metal_node
        row["organic_linker"] = organic_linker
        normalized_rows.append({field: row.get(field, "") for field in FIELDS})
    rows = normalized_rows
    rows.sort(key=lambda row: (row["basename"], int(row["temperature_K"]), int(row["pressure_bar"])))
    with FINAL.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {FINAL}: rows={len(rows)}, completed={sum(r['completed'] == 'yes' for r in rows)}, parsed={sum(r['result_parse_status'] == 'parsed' for r in rows)}")


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--temperature", type=int)
    group.add_argument("--combine", action="store_true")
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=4)
    args = parser.parse_args()
    if args.shard_count < 1 or not 0 <= args.shard_index < args.shard_count:
        parser.error("shard index/count are inconsistent")
    combine(args.shard_count) if args.combine else process_temperature(args.temperature, args.shard_index, args.shard_count)


if __name__ == "__main__":
    main()
