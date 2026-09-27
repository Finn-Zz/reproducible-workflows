#!/usr/bin/env python3
"""Rebuild v5 from a verified second-round archive and a v4 summary.

This is a provenance/rebuild utility; the already validated v5 CSV is shipped
under ``data/`` and is the normal analysis input.  The raw ZIP and v4 table are
read-only. A job is accepted only when its current top-level output contains
``Simulation finished!`` and all final loading fields can be parsed. Status TSVs
and success markers are cross-checks, not substitutes for a parseable output.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import hashlib
import io
import math
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


# ``base`` is loaded from an explicit path in ``main``.  Keeping the raw
# archive and v4 table outside this small package makes the package portable.
base = None
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
WORK_ROOT = PACKAGE_ROOT
DEFAULT_RESULT_ROOT = PACKAGE_ROOT / "raw"
DEFAULT_BASE_SCRIPT = PACKAGE_ROOT / "scripts" / "build_kisti_overall_summary.py"
ARCHIVE = DEFAULT_RESULT_ROOT / "2round (1).zip"
SOURCE = PACKAGE_ROOT / "data" / "KISTI_adsorption_overall_summary_v4.csv"
TARGET = PACKAGE_ROOT / "data" / "KISTI_adsorption_overall_summary_v5_rebuilt.csv"
REPORT = PACKAGE_ROOT / "data" / "2round_1_update_v5_rebuilt_report.md"


def load_base_module(path: Path):
    """Load the parser used for the original summary without PYTHONPATH hacks."""
    if not path.is_file():
        raise FileNotFoundError(
            f"Parser module not found: {path}. Supply --base-script or use the bundled copy."
        )
    spec = importlib.util.spec_from_file_location("kisti_summary_parser", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load parser module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise RuntimeError(f"Missing CSV header: {path}")
        return reader.fieldnames, list(reader)


def write_csv_atomic(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def table_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return row["temperature_K"], row["topology"], row["pressure_bar"], row["basename"]


def status_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return row["temperature_K"], row["topology"], row["pressure_bar"], row["basename"]


def output_candidates(status: dict[str, str]) -> tuple[str, ...]:
    temperature = status["temperature_K"]
    pressure_token = "5e+05" if status["pressure_bar"] == "5" else "1e+07"
    filename = f"output_{temperature}_{pressure_token}.s0.txt"
    prefix = status["relative_job_directory"].rstrip("/")
    return (
        f"{prefix}/output_{temperature}/{filename}",
        f"{prefix}/output/{filename}",
        f"{prefix}/{filename}",
    )


def simulation_member(status: dict[str, str]) -> str:
    return status["relative_job_directory"].rstrip("/") + "/simulation.json"


def validate_final_output(
    text: str,
    parsed: dict[str, str],
    temperature: int,
    pressure_bar: int,
    member: str,
) -> None:
    missing = [field for field in base.LOADING_FIELDS if field not in parsed]
    if "Simulation finished!" not in text or missing:
        raise RuntimeError(f"Incomplete final output classified as complete: {member}; missing={missing}")
    for field, expected in (("Temperature", temperature), ("Pressure", pressure_bar * 100000)):
        match = re.search(rf"^{field}:\s*({base.NUMBER})", text, re.MULTILINE)
        if not match or not math.isclose(float(match.group(1)), expected, rel_tol=0.0, abs_tol=1e-7):
            raise RuntimeError(f"{member}: incorrect or missing {field}")
    if not math.isclose(
        float(parsed["absolute_loading_mg_g_framework"]),
        float(parsed["absolute_loading_mol_kg_framework"]) * 2.016,
        rel_tol=2e-5,
    ):
        raise RuntimeError(f"{member}: inconsistent absolute-loading units")
    if float(parsed["absolute_loading_mg_g_framework"]) < 0:
        raise RuntimeError(f"{member}: negative absolute loading")
    if not all(math.isfinite(float(value)) for value in parsed.values()):
        raise RuntimeError(f"{member}: non-finite parsed value")


def main() -> int:
    global ARCHIVE, SOURCE, TARGET, REPORT, base
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ARCHIVE,
                        help=f"verified second-round ZIP (default: {ARCHIVE})")
    parser.add_argument("--source", type=Path, default=SOURCE,
                        help=f"v4 summary table (default: {SOURCE})")
    parser.add_argument("--target", type=Path, default=TARGET,
                        help=f"new v5 output table (default: {TARGET})")
    parser.add_argument("--report", type=Path, default=REPORT,
                        help=f"merge report (default: {REPORT})")
    parser.add_argument("--base-script", type=Path, default=DEFAULT_BASE_SCRIPT,
                        help=f"summary parser module (default: {DEFAULT_BASE_SCRIPT})")
    parser.add_argument("--overwrite", action="store_true", help="allow replacing --target and --report")
    args = parser.parse_args()
    ARCHIVE, SOURCE, TARGET, REPORT = (
        args.archive.expanduser().resolve(), args.source.expanduser().resolve(),
        args.target.expanduser().resolve(), args.report.expanduser().resolve(),
    )
    base = load_base_module(args.base_script.expanduser().resolve())
    if TARGET.exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite existing output: {TARGET}; use --overwrite")
    if not ARCHIVE.is_file() or not SOURCE.is_file():
        raise FileNotFoundError(
            f"Required ZIP or v4 source table is missing: archive={ARCHIVE}, source={SOURCE}"
        )
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    REPORT.parent.mkdir(parents=True, exist_ok=True)

    source_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    archive_hash = hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()
    fields, source_rows = read_csv(SOURCE)
    if fields != ["current_basename", *base.FIELDS]:
        raise RuntimeError(f"Unexpected v4 schema: {fields}")
    indexed = {table_key(row): row.copy() for row in source_rows}
    if len(source_rows) != 6912 or len(indexed) != len(source_rows):
        raise RuntimeError(f"Unexpected v4 row/key count: {len(source_rows)}/{len(indexed)}")

    updated_keys: list[tuple[str, str, str, str]] = []
    verified_existing_keys: list[tuple[str, str, str, str]] = []
    archive_incomplete: list[dict[str, str]] = []
    per_temperature: dict[int, Counter[str]] = defaultdict(Counter)

    with zipfile.ZipFile(ARCHIVE) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise RuntimeError("ZIP contains duplicate member names")
        status_members = sorted(
            name for name in names
            if re.fullmatch(r"check_continue_status_\d{8}T\d{6}\.tsv", name)
            and not name.endswith("205021.tsv")
        )
        # The final pre-submission snapshot is authoritative for the 153-job registry.
        completed_snapshots = [name for name in status_members if name.endswith("204833.tsv")]
        if len(completed_snapshots) != 1:
            raise RuntimeError(f"Expected one final diagnostic status snapshot, found {completed_snapshots}")
        status_member = completed_snapshots[0]
        statuses = list(csv.DictReader(
            io.StringIO(archive.read(status_member).decode("utf-8")), delimiter="\t"
        ))
        if len(statuses) != 153 or len({status_key(row) for row in statuses}) != 153:
            raise RuntimeError("Final status snapshot does not contain 153 unique jobs")

        actual_complete = 0
        for status in statuses:
            key = status_key(status)
            if key not in indexed:
                raise KeyError(f"2round task is absent from v4: {key}")
            sim_member = simulation_member(status)
            if sim_member not in archive.NameToInfo:
                raise RuntimeError(f"Missing simulation.json: {sim_member}")
            output_member = next(
                (candidate for candidate in output_candidates(status) if candidate in archive.NameToInfo), None
            )
            text = archive.read(output_member).decode("utf-8", errors="replace") if output_member else ""
            parsed = base.parse_output(text) if text else {}
            complete = bool(
                output_member
                and "Simulation finished!" in text
                and all(field in parsed for field in base.LOADING_FIELDS)
            )
            status_says_complete = status["status"] == "complete"
            success_marker = (
                status["relative_job_directory"].rstrip("/")
                + f"/.raspa_success_{status['temperature_K']}"
            ) in archive.NameToInfo
            if complete != status_says_complete or complete != success_marker:
                raise RuntimeError(
                    f"Completion evidence disagrees for {key}: output={complete}, "
                    f"status={status['status']}, marker={success_marker}"
                )
            per_temperature[int(status["temperature_K"])]["complete" if complete else "incomplete"] += 1
            if not complete:
                archive_incomplete.append(status)
                continue

            actual_complete += 1
            validate_final_output(
                text, parsed, int(status["temperature_K"]), int(status["pressure_bar"]), output_member
            )
            row = indexed[key]
            if row["completed"] == "yes":
                # These 59 results were already applied to v4.  Confirm numerical identity and leave them bytewise unchanged.
                for field in (*base.LOADING_FIELDS, *base.EXTRA_FIELDS):
                    if row.get(field, "") and parsed.get(field, "") and not math.isclose(
                        float(row[field]), float(parsed[field]), rel_tol=1e-12, abs_tol=0.0
                    ):
                        raise RuntimeError(f"Previously merged value changed for {key}, field {field}")
                verified_existing_keys.append(key)
            else:
                row.update(parsed)
                row["completed"] = "yes"
                row["completed_via"] = "second"
                row["result_parse_status"] = "parsed"
                updated_keys.append(key)

        if actual_complete != 126 or len(verified_existing_keys) != 59 or len(updated_keys) != 67:
            raise RuntimeError(
                f"Unexpected merge counts: complete={actual_complete}, existing={len(verified_existing_keys)}, "
                f"new={len(updated_keys)}"
            )
        if len(archive_incomplete) != 27:
            raise RuntimeError(f"Unexpected incomplete ZIP-job count: {len(archive_incomplete)}")

    target_rows = [indexed[table_key(row)] for row in source_rows]
    if hashlib.sha256(SOURCE.read_bytes()).hexdigest() != source_hash:
        raise RuntimeError("v4 changed while v5 was being generated")
    untouched = set(indexed) - set(updated_keys)
    source_index = {table_key(row): row for row in source_rows}
    if any(indexed[key] != source_index[key] for key in untouched):
        raise RuntimeError("A v4 row outside the 67 intended updates was changed")
    complete_total = sum(row["completed"] == "yes" for row in target_rows)
    incomplete_rows = [row for row in target_rows if row["completed"] != "yes"]
    if complete_total != 6822 or len(incomplete_rows) != 90:
        raise RuntimeError(f"Unexpected v5 completion totals: {complete_total} complete, {len(incomplete_rows)} incomplete")
    if len({row["current_basename"] for row in target_rows}) != 432:
        raise RuntimeError("Unexpected v5 structure count")

    write_csv_atomic(TARGET, fields, target_rows)
    target_fields, roundtrip = read_csv(TARGET)
    if target_fields != fields or roundtrip != target_rows:
        raise RuntimeError("v5 round-trip verification failed")

    overall_by_temp: dict[int, Counter[str]] = defaultdict(Counter)
    remaining_by_status = Counter()
    remaining_by_temperature = Counter()
    for row in target_rows:
        state = "complete" if row["completed"] == "yes" else "incomplete"
        overall_by_temp[int(row["temperature_K"])][state] += 1
        if state == "incomplete":
            remaining_by_status[row["result_parse_status"]] += 1
            remaining_by_temperature[int(row["temperature_K"])] += 1

    report_lines = [
        "# `2round (1).zip` 合并到 v5 的核对报告",
        "",
        "- 基础表：`KISTI_adsorption_overall_summary_v4.csv`（未覆盖）。",
        "- 新表：`KISTI_adsorption_overall_summary_v5.csv`。",
        "- ZIP登记153个二轮任务：126个具有完成标志和可完整解析的终态输出，27个仍未完成。",
        "- 126个完成任务中，59个已经存在于v4且数值一致；此次新增67个完成结果。",
        "- v5共6912行、432个结构：6822个任务完成，90个任务未完成。",
        "- v5字段与v4完全一致，没有增加误差列。",
        "",
        "## v5总体完成情况",
        "",
        "| 温度/K | 完成 | 未完成 |",
        "|---:|---:|---:|",
    ]
    for temperature in sorted(overall_by_temp):
        counts = overall_by_temp[temperature]
        report_lines.append(f"| {temperature} | {counts['complete']} | {counts['incomplete']} |")
    report_lines += [
        "",
        "## 本ZIP内任务",
        "",
        "| 温度/K | ZIP内完成 | ZIP内未完成 |",
        "|---:|---:|---:|",
    ]
    for temperature in sorted(per_temperature):
        counts = per_temperature[temperature]
        report_lines.append(f"| {temperature} | {counts['complete']} | {counts['incomplete']} |")
    report_lines += [
        "",
        "## v5尚未完成任务的原状态",
        "",
        "| result_parse_status | 数量 |",
        "|---|---:|",
    ]
    for status, count in sorted(remaining_by_status.items()):
        report_lines.append(f"| {status} | {count} |")
    report_lines += [
        "",
        "按温度统计尚未完成："
        + "；".join(f"{temperature} K：{count}" for temperature, count in sorted(remaining_by_temperature.items()))
        + "。",
        "",
        "完成判据：当前顶层输出同时包含 `Simulation finished!` 和全部最终吸附量字段；单独的状态表或success标记不作为完成依据。",
        "",
        f"v4 SHA256：`{source_hash}`",
        f"ZIP SHA256：`{archive_hash}`",
        f"v5 SHA256：`{hashlib.sha256(TARGET.read_bytes()).hexdigest()}`",
        "",
    ]
    REPORT.write_text("\n".join(report_lines), encoding="utf-8")

    print(f"Verified archive jobs: 153 (126 complete, 27 incomplete)")
    print(f"Newly applied to v4: {len(updated_keys)}")
    print(f"Already present and numerically identical: {len(verified_existing_keys)}")
    print(f"v5: {len(target_rows)} rows, {complete_total} complete, {len(incomplete_rows)} incomplete")
    print(f"Wrote: {TARGET}")
    print(f"Report: {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
