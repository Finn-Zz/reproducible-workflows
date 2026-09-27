#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(os.environ.get("HMOF_WORKSPACE_ROOT", "."))
STRUCTURES = ROOT / "Structures_MS"
LI_PYTHON = Path(os.environ.get("LAMMPS_INTERFACE_PYTHON", sys.executable))
METALS = ("F", "TiCo", "TiMg", "TiNi")

ATOM_MINIMIZE = "minimize        1.0e-6 1.0e-6 20000 200000"
CELL_MINIMIZE_BLOCK = [
    "fix             cellrelax all box/relax aniso 0.0 vmax 0.0005",
    "minimize        1.0e-6 1.0e-6 20000 200000",
    "unfix           cellrelax",
    "minimize        1.0e-6 1.0e-6 10000 100000",
]


def run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        cmd,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def rename_ms_dirs() -> None:
    for metal in METALS:
        old = STRUCTURES / f"{metal}_2"
        new = STRUCTURES / f"{metal}_2_MS"
        if new.exists():
            continue
        if old.exists():
            old.rename(new)


def convert_cif(cif: Path, workdir: Path) -> subprocess.CompletedProcess[str]:
    code = (
        "import sys\n"
        "sys.setrecursionlimit(50000)\n"
        "from lammps_interface.cli import main\n"
        "raise SystemExit(main())\n"
    )
    return run(
        [str(LI_PYTHON), "-c", code, str(cif), "-ff", "UFF4MOF", "--relax"],
        cwd=workdir,
    )


def patch_common(line: str, name: str, suffix: str) -> str | None:
    stripped = line.strip()
    if stripped.startswith("print "):
        return None
    if stripped.startswith("log "):
        return f"log             log.{name}.{suffix} append"
    if stripped.startswith("dump "):
        parts = line.split()
        if len(parts) >= 7:
            parts[1] = f"{name}_{suffix}"
            parts[5] = f"relaxed_{name}_{suffix}.lammpstrj"
            return " ".join(parts)
    if stripped.startswith("dump_modify "):
        return line.replace(f"{name}_relax", f"{name}_{suffix}", 1)
    return line


def make_atom_input(raw_input: Path, out_input: Path, name: str) -> None:
    lines: list[str] = []
    for line in raw_input.read_text().splitlines():
        patched = patch_common(line, name, "atom")
        if patched is None:
            continue
        if patched.strip().startswith("minimize "):
            patched = ATOM_MINIMIZE
        lines.append(patched)
    out_input.write_text("\n".join(lines) + "\n")


def make_cell_input(raw_input: Path, out_input: Path, name: str) -> None:
    lines: list[str] = []
    inserted = False
    for line in raw_input.read_text().splitlines():
        patched = patch_common(line, name, "cell")
        if patched is None:
            continue
        if patched.strip().startswith("minimize "):
            lines.extend(CELL_MINIMIZE_BLOCK)
            inserted = True
            continue
        lines.append(patched)
    if not inserted:
        raise RuntimeError(f"No minimize line found in {raw_input}")
    out_input.write_text("\n".join(lines) + "\n")


def prepare_one(cif: Path, out_root: Path) -> dict[str, str]:
    name = cif.stem
    workdir = out_root / name
    workdir.mkdir(parents=True, exist_ok=True)
    copied_cif = workdir / cif.name
    if not copied_cif.exists():
        shutil.copy2(cif, copied_cif)

    row = {
        "name": name,
        "source_cif": str(cif.relative_to(ROOT)),
        "workdir": str(workdir.relative_to(ROOT)),
        "conversion_status": "",
        "atom_input": "",
        "cell_input": "",
        "error_message": "",
    }

    raw_input = workdir / f"in.{name}.raw"
    atom_input = workdir / f"in.{name}.atom"
    cell_input = workdir / f"in.{name}.cell"
    data_file = workdir / f"data.{name}"
    if atom_input.exists() and cell_input.exists() and data_file.exists():
        row.update(
            {
                "conversion_status": "ok_existing",
                "atom_input": atom_input.name,
                "cell_input": cell_input.name,
            }
        )
        return row

    conv = convert_cif(copied_cif, workdir)
    (workdir / "lammps_interface.stdout").write_text(conv.stdout)
    if conv.returncode != 0:
        row["conversion_status"] = f"failed:{conv.returncode}"
        row["error_message"] = conv.stdout.splitlines()[-1] if conv.stdout.splitlines() else ""
        return row

    generated = workdir / f"in.{name}"
    if not generated.exists():
        row["conversion_status"] = "failed:missing_in"
        row["error_message"] = f"Missing generated {generated.name}"
        return row
    generated.replace(raw_input)

    try:
        make_atom_input(raw_input, atom_input, name)
        make_cell_input(raw_input, cell_input, name)
    except Exception as exc:
        row["conversion_status"] = "failed:patch"
        row["error_message"] = repr(exc)
        return row

    row.update(
        {
            "conversion_status": "ok",
            "atom_input": atom_input.name,
            "cell_input": cell_input.name,
        }
    )
    return row


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    rename_ms_dirs()
    all_rows: list[dict[str, str]] = []
    for metal in METALS:
        source = STRUCTURES / f"{metal}_2_MS"
        out_root = STRUCTURES / f"{metal}_3_lammps"
        out_root.mkdir(parents=True, exist_ok=True)
        cifs = sorted(source.glob("*.cif"))
        rows: list[dict[str, str]] = []
        for idx, cif in enumerate(cifs, 1):
            print(f"[{metal} {idx}/{len(cifs)}] {cif.name}", flush=True)
            row = prepare_one(cif, out_root)
            rows.append(row)
            all_rows.append({"metal": metal, **row})
            write_summary(out_root / f"{metal}_3_lammps_conversion_summary.csv", rows)
            write_summary(STRUCTURES / "lammps_conversion_summary.csv", all_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
