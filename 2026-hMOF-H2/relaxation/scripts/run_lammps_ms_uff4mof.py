#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import re
import subprocess
from pathlib import Path

from ase.io import read, write


ROOT = Path(os.environ.get("HMOF_WORKSPACE_ROOT", "."))
STRUCTURES = ROOT / "Structures_MS"
METALS = ("F", "TiCo", "TiMg", "TiNi")
STEPS = ("atom", "cell")
LMP = os.environ.get("LMP_BIN", "lmp")
THREADS = os.environ.get("LAMMPS_OMP_THREADS", "8")
SUMMARY = STRUCTURES / "lammps_ms_uff4mof_optimization_summary.csv"


def parse_lammps_output(text: str) -> dict[str, str]:
    errors = re.findall(r"ERROR[^\n]*(?:\n[^\n]*)?", text)
    forces = re.findall(
        r"Force max component initial, final =\s+([\deE+\-.]+)\s+([\deE+\-.]+)",
        text,
    )
    criteria = re.findall(r"Stopping criterion =\s+([^\n]+)", text)
    energies = re.findall(
        r"Energy initial, next-to-last, final =\s+([\deE+\-.]+)\s+([\deE+\-.]+)\s+([\deE+\-.]+)",
        text,
    )
    return {
        "lammps_error": "yes" if errors else "no",
        "error_message": errors[-1].replace("\n", " | ") if errors else "",
        "initial_force_max": forces[-1][0] if forces else "",
        "final_force_max": forces[-1][1] if forces else "",
        "initial_energy": energies[-1][0] if energies else "",
        "final_energy": energies[-1][2] if energies else "",
        "stopping_criterion": criteria[-1].strip() if criteria else "",
    }


def write_rows(rows: list[dict[str, str]]) -> None:
    fields = [
        "metal",
        "name",
        "step",
        "workdir",
        "input_file",
        "data_file",
        "lammps_status",
        "has_output_cif",
        "output_cif",
        "lammps_error",
        "initial_force_max",
        "final_force_max",
        "initial_energy",
        "final_energy",
        "stopping_criterion",
        "error_message",
    ]
    with SUMMARY.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def existing_rows() -> list[dict[str, str]]:
    if not SUMMARY.exists():
        return []
    with SUMMARY.open(newline="") as fh:
        return list(csv.DictReader(fh))


def row_key(row: dict[str, str]) -> tuple[str, str, str]:
    return row["metal"], row["name"], row["step"]


def run_one(metal: str, workdir: Path, step: str) -> dict[str, str]:
    name = workdir.name
    input_file = workdir / f"in.{name}.{step}"
    data_file = workdir / f"data.{name}"
    stdout_file = workdir / f"lammps_{step}.stdout"
    dump_file = workdir / f"relaxed_{name}_{step}.lammpstrj"
    output_cif = workdir / f"{name}_UFF4MOF_{step}.cif"

    row = {
        "metal": metal,
        "name": name,
        "step": step,
        "workdir": str(workdir.relative_to(ROOT)),
        "input_file": input_file.name,
        "data_file": data_file.name,
        "lammps_status": "",
        "has_output_cif": "no",
        "output_cif": "",
        "lammps_error": "",
        "initial_force_max": "",
        "final_force_max": "",
        "initial_energy": "",
        "final_energy": "",
        "stopping_criterion": "",
        "error_message": "",
    }

    if not input_file.exists() or not data_file.exists():
        row["lammps_status"] = "missing_input"
        row["error_message"] = "missing LAMMPS input or data file"
        return row

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = THREADS
    cmd = [LMP, "-sf", "omp", "-pk", "omp", THREADS, "-in", input_file.name]
    proc = subprocess.run(
        cmd,
        cwd=workdir,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    stdout_file.write_text(proc.stdout)
    parsed = parse_lammps_output(proc.stdout)
    row.update(parsed)
    row["lammps_status"] = "ok" if proc.returncode == 0 and parsed["lammps_error"] == "no" else f"failed:{proc.returncode}"

    if dump_file.exists() and row["lammps_status"] == "ok":
        try:
            write(output_cif, read(dump_file, index=-1))
            row["has_output_cif"] = "yes"
            row["output_cif"] = output_cif.name
        except Exception as exc:
            row["lammps_status"] = "failed:ase_convert"
            row["error_message"] = str(exc)

    return row


def main() -> int:
    rows = existing_rows()
    by_key = {row_key(row): row for row in rows}
    tasks: list[tuple[str, Path, str]] = []
    for metal in METALS:
        root = STRUCTURES / f"{metal}_3_lammps"
        for workdir in sorted(p for p in root.iterdir() if p.is_dir()):
            for step in STEPS:
                key = (metal, workdir.name, step)
                row = by_key.get(key)
                output_cif = workdir / f"{workdir.name}_UFF4MOF_{step}.cif"
                if row and row.get("lammps_status") == "ok" and output_cif.exists():
                    continue
                tasks.append((metal, workdir, step))

    total = len(tasks)
    print(f"pending tasks: {total}", flush=True)
    for idx, (metal, workdir, step) in enumerate(tasks, 1):
        print(f"[{idx}/{total}] {metal} {step} {workdir.name}", flush=True)
        row = run_one(metal, workdir, step)
        by_key[(metal, workdir.name, step)] = row
        rows = [by_key[key] for key in sorted(by_key)]
        write_rows(rows)
        print(
            f"    {row['lammps_status']} force={row['final_force_max']} criterion={row['stopping_criterion']}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
