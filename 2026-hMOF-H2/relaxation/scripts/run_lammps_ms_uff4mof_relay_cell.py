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
LMP = os.environ.get("LMP_BIN", "lmp")
THREADS = os.environ.get("LAMMPS_OMP_THREADS", "8")
SUMMARY = STRUCTURES / "lammps_ms_uff4mof_relay_atom_cell_summary.csv"


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


def run_lmp(workdir: Path, input_name: str) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = THREADS
    return subprocess.run(
        [LMP, "-sf", "omp", "-pk", "omp", THREADS, "-in", input_name],
        cwd=workdir,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def clean_old_cell_outputs(workdir: Path, name: str) -> None:
    patterns = [
        f"{name}_UFF4MOF_cell.cif",
        f"{name}_UFF4MOF_cell_from_atom.cif",
        f"relaxed_{name}_cell.lammpstrj",
        f"relaxed_{name}_cell_from_atom.lammpstrj",
        "lammps_cell.stdout",
        "lammps_cell_from_atom.stdout",
        f"log.{name}.cell",
        f"log.{name}.cell_from_atom",
        f"data.{name}.atom_relaxed",
        f"restart.{name}.atom_relaxed",
        f"in.{name}.atom_to_data",
        f"in.{name}.atom_to_restart",
        f"in.{name}.atom_to_dump",
        f"in.{name}.cell_from_atom",
        f"atom_relaxed_{name}_for_cell.lammpstrj",
        "lammps_atom_to_dump.stdout",
    ]
    for pattern in patterns:
        for path in workdir.glob(pattern):
            if path.is_file():
                path.unlink()


CELL_MINIMIZE_BLOCK = [
    "fix             cellrelax all box/relax aniso 0.0 vmax 0.0005",
    "minimize        1.0e-6 1.0e-6 20000 200000",
    "unfix           cellrelax",
    "minimize        1.0e-6 1.0e-6 10000 100000",
]


def make_atom_to_dump_input(workdir: Path, name: str) -> Path:
    src = workdir / f"in.{name}.atom"
    dst = workdir / f"in.{name}.atom_to_dump"
    lines: list[str] = []
    for line in src.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("log "):
            lines.append(f"log             log.{name}.atom_to_dump append")
            continue
        if stripped.startswith("dump ") or stripped.startswith("dump_modify "):
            continue
        lines.append(line)
    lines.append("reset_timestep  0")
    lines.append(
        f"write_dump      all custom atom_relaxed_{name}_for_cell.lammpstrj "
        "id type x y z modify sort id"
    )
    dst.write_text("\n".join(lines) + "\n")
    return dst


def make_cell_from_atom_input(workdir: Path, name: str) -> Path:
    src = workdir / f"in.{name}.raw"
    dst = workdir / f"in.{name}.cell_from_atom"
    lines: list[str] = []
    inserted_cell = False
    for line in src.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("print "):
            continue
        if stripped.startswith("log "):
            lines.append(f"log             log.{name}.cell_from_atom append")
            continue
        if stripped.startswith("read_data "):
            lines.append(f"read_data       data.{name}")
            lines.append(
                f"read_dump       atom_relaxed_{name}_for_cell.lammpstrj 0 "
                "x y z box yes replace yes"
            )
            continue
        if stripped.startswith("dump "):
            parts = line.split()
            if len(parts) >= 7:
                parts[1] = f"{name}_cell_from_atom"
                parts[5] = f"relaxed_{name}_cell_from_atom.lammpstrj"
                lines.append(" ".join(parts))
                continue
        if stripped.startswith("dump_modify "):
            parts = line.split()
            if len(parts) >= 2:
                parts[1] = f"{name}_cell_from_atom"
                lines.append(" ".join(parts))
                continue
        if stripped.startswith("minimize "):
            lines.extend(CELL_MINIMIZE_BLOCK)
            inserted_cell = True
            continue
        lines.append(line)
    if not inserted_cell:
        raise RuntimeError(f"No minimize line found in {src}")
    dst.write_text("\n".join(lines) + "\n")
    return dst


def write_rows(rows: list[dict[str, str]]) -> None:
    fields = [
        "metal",
        "name",
        "stage",
        "workdir",
        "input_file",
        "output_data",
        "output_cif",
        "lammps_status",
        "has_output",
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


def existing_rows() -> dict[tuple[str, str, str], dict[str, str]]:
    if not SUMMARY.exists():
        return {}
    with SUMMARY.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    return {(r["metal"], r["name"], r["stage"]): r for r in rows}


def atom_to_dump(metal: str, workdir: Path) -> dict[str, str]:
    name = workdir.name
    input_file = make_atom_to_dump_input(workdir, name)
    output_dump = workdir / f"atom_relaxed_{name}_for_cell.lammpstrj"
    proc = run_lmp(workdir, input_file.name)
    (workdir / "lammps_atom_to_dump.stdout").write_text(proc.stdout)
    parsed = parse_lammps_output(proc.stdout)
    status = "ok" if proc.returncode == 0 and parsed["lammps_error"] == "no" and output_dump.exists() else f"failed:{proc.returncode}"
    return {
        "metal": metal,
        "name": name,
        "stage": "atom_to_dump",
        "workdir": str(workdir.relative_to(ROOT)),
        "input_file": input_file.name,
        "output_data": output_dump.name if output_dump.exists() else "",
        "output_cif": "",
        "lammps_status": status,
        "has_output": "yes" if output_dump.exists() else "no",
        **parsed,
    }


def cell_from_atom(metal: str, workdir: Path) -> dict[str, str]:
    name = workdir.name
    input_file = make_cell_from_atom_input(workdir, name)
    dump = workdir / f"relaxed_{name}_cell_from_atom.lammpstrj"
    output_cif = workdir / f"{name}_UFF4MOF_cell.cif"
    proc = run_lmp(workdir, input_file.name)
    (workdir / "lammps_cell_from_atom.stdout").write_text(proc.stdout)
    parsed = parse_lammps_output(proc.stdout)
    status = "ok" if proc.returncode == 0 and parsed["lammps_error"] == "no" else f"failed:{proc.returncode}"
    has_output = "no"
    if dump.exists() and status == "ok":
        try:
            write(output_cif, read(dump, index=-1))
            has_output = "yes"
        except Exception as exc:
            status = "failed:ase_convert"
            parsed["error_message"] = str(exc)
    return {
        "metal": metal,
        "name": name,
        "stage": "cell_from_atom",
        "workdir": str(workdir.relative_to(ROOT)),
        "input_file": input_file.name,
        "output_data": "",
        "output_cif": output_cif.name if output_cif.exists() else "",
        "lammps_status": status,
        "has_output": has_output,
        **parsed,
    }


def main() -> int:
    by_key = existing_rows()
    tasks = [(metal, d) for metal in METALS for d in sorted((STRUCTURES / f"{metal}_3_lammps").iterdir()) if d.is_dir()]
    total = len(tasks)
    for idx, (metal, workdir) in enumerate(tasks, 1):
        name = workdir.name
        atom_key = (metal, name, "atom_to_dump")
        cell_key = (metal, name, "cell_from_atom")
        print(f"[{idx}/{total}] {metal} {name}", flush=True)
        if atom_key not in by_key or by_key[atom_key].get("lammps_status") != "ok":
            clean_old_cell_outputs(workdir, name)
            row = atom_to_dump(metal, workdir)
            by_key[atom_key] = row
            write_rows([by_key[k] for k in sorted(by_key)])
            print(f"    atom_to_dump {row['lammps_status']} force={row['final_force_max']}", flush=True)
            if row["lammps_status"] != "ok":
                continue
        if cell_key not in by_key or by_key[cell_key].get("lammps_status") != "ok":
            row = cell_from_atom(metal, workdir)
            by_key[cell_key] = row
            write_rows([by_key[k] for k in sorted(by_key)])
            print(f"    cell_from_atom {row['lammps_status']} force={row['final_force_max']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
