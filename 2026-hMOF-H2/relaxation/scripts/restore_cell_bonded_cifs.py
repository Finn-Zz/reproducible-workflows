#!/usr/bin/env python3
"""Restore LAMMPS cell-relaxed topology and write bonded P1 CIF files.

ASE CIF export retains coordinates and cell but not the LAMMPS Bonds section.
This script transfers the final scaled coordinates from the cell-relax dump onto
the original LAMMPS data topology.  It never infers bonds from distances.
"""

from __future__ import annotations

import argparse
import math
import os
import re
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment


ROOT = Path(os.environ.get("HMOF_STRUCTURES_ROOT", Path(os.environ.get("HMOF_WORKSPACE_ROOT", ".")) / "Structures_MS"))
METALS = ("F", "TiCo", "TiMg", "TiNi")
ELEMENT_RE = re.compile(r"^([A-Z][a-z]?)")


def parse_data(path: Path):
    lines = path.read_text().splitlines()
    atom_types: dict[int, str] = {}
    bonds: list[tuple[int, int, int]] = []
    atoms_start = bonds_start = None
    for i, line in enumerate(lines):
        if line.strip() == "Masses":
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            while j < len(lines) and lines[j].strip():
                fields = lines[j].split("#", 1)
                atom_type = int(fields[0].split()[0])
                label = fields[1].strip() if len(fields) == 2 else "X"
                match = ELEMENT_RE.match(label)
                atom_types[atom_type] = match.group(1) if match else "X"
                j += 1
        if line.strip() == "Atoms":
            atoms_start = i + 1
        if line.strip() == "Bonds":
            bonds_start = i + 1

    if atoms_start is None or bonds_start is None:
        raise ValueError(f"Missing Atoms or Bonds section: {path}")

    atoms: dict[int, tuple[int, int]] = {}
    i = atoms_start
    while i < len(lines) and not lines[i].strip():
        i += 1
    while i < len(lines) and lines[i].strip():
        values = lines[i].split()
        atoms[int(values[0])] = (i, int(values[2]))
        i += 1

    i = bonds_start
    while i < len(lines) and not lines[i].strip():
        i += 1
    while i < len(lines) and lines[i].strip():
        values = lines[i].split()
        bonds.append((int(values[1]), int(values[2]), int(values[3])))
        i += 1
    return lines, atoms, atom_types, bonds


def parse_box(bounds):
    """Return the restricted triclinic cell and true LAMMPS box origin."""
    xlo_bound, xhi_bound, xy = bounds[0]
    ylo_bound, yhi_bound, xz = bounds[1]
    zlo_bound, zhi_bound, yz = bounds[2]
    xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
    xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
    ylo = ylo_bound - min(0.0, yz)
    yhi = yhi_bound - max(0.0, yz)
    cell = (
        (xhi - xlo, 0.0, 0.0),
        (xy, yhi - ylo, 0.0),
        (xz, yz, zhi_bound - zlo_bound),
    )
    return (xlo, ylo, zlo_bound), cell


def final_frame(path: Path):
    """Read the final LAMMPS custom dump frame."""
    lines = path.read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if line == "ITEM: TIMESTEP"]
    if not starts:
        raise ValueError(f"No dump frame found: {path}")
    i = starts[-1]
    n_atoms = int(lines[i + 3])
    if not lines[i + 4].startswith("ITEM: BOX BOUNDS"):
        raise ValueError(f"Unexpected box block: {path}")
    bounds = [[float(x) for x in lines[i + 5 + j].split()] for j in range(3)]
    atoms_header = lines[i + 8].split()[2:]
    if atoms_header != ["element", "xs", "ys", "zs"]:
        raise ValueError(f"Unexpected dump columns {atoms_header}: {path}")
    coordinates = []
    elements = []
    for line in lines[i + 9 : i + 9 + n_atoms]:
        element, xs, ys, zs = line.split()
        elements.append(element)
        coordinates.append((float(xs), float(ys), float(zs)))
    if len(coordinates) != n_atoms:
        raise ValueError(f"Truncated final frame: {path}")

    origin, cell = parse_box(bounds)
    return elements, coordinates, cell, origin


def parse_atom_relaxed_dump(path: Path, atom_types):
    """Read the ID-labelled atom-relaxed structure used as the cell input."""
    lines = path.read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if line == "ITEM: TIMESTEP"]
    if not starts:
        raise ValueError(f"No dump frame found: {path}")
    i = starts[-1]
    n_atoms = int(lines[i + 3])
    bounds = [[float(x) for x in lines[i + 5 + j].split()] for j in range(3)]
    header = lines[i + 8].split()[2:]
    if header != ["id", "type", "x", "y", "z"]:
        raise ValueError(f"Unexpected atom-relaxed columns {header}: {path}")
    origin, cell = parse_box(bounds)
    inverse_cell_transpose = np.linalg.inv(np.asarray(cell).T)
    frac_by_id = {}
    element_by_id = {}
    for line in lines[i + 9 : i + 9 + n_atoms]:
        atom_id, atom_type, x, y, z = line.split()
        atom_id = int(atom_id)
        position = np.asarray((float(x), float(y), float(z))) - np.asarray(origin)
        frac_by_id[atom_id] = tuple((inverse_cell_transpose @ position) % 1.0)
        element_by_id[atom_id] = atom_types[int(atom_type)]
    if len(frac_by_id) != n_atoms:
        raise ValueError(f"Truncated atom-relaxed frame: {path}")
    return frac_by_id, element_by_id


def cartesian(frac, cell):
    return tuple(sum(frac[j] * cell[j][i] for j in range(3)) for i in range(3))


def cell_parameters(cell):
    def norm(v):
        return math.sqrt(sum(x * x for x in v))

    def angle(u, v):
        cosine = sum(a * b for a, b in zip(u, v)) / (norm(u) * norm(v))
        return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))

    a, b, c = cell
    return norm(a), norm(b), norm(c), angle(b, c), angle(a, c), angle(a, b)


def nearest_image_delta(frac1, frac2, cell):
    best = None
    for tx in (-1, 0, 1):
        for ty in (-1, 0, 1):
            for tz in (-1, 0, 1):
                delta = tuple(frac2[k] + (tx, ty, tz)[k] - frac1[k] for k in range(3))
                vector = cartesian(delta, cell)
                distance = math.sqrt(sum(x * x for x in vector))
                candidate = (distance, (tx, ty, tz))
                if best is None or candidate < best:
                    best = candidate
    assert best is not None
    return best


def periodic_code(translation):
    # CIF symmetry code: identity operation 1 plus translation digits offset by 5.
    return "1_" + "".join(str(component + 5) for component in translation)


def match_final_coordinates(start_coords, start_elements, final_elements, final_coords, final_cell):
    """Map the unlabeled final dump back to atom IDs by minimum PBC displacement.

    LAMMPS reorders atoms internally during minimization, while the final dump
    written by the historical workflow did not include `id`.  The preceding
    atom-relaxed dump does include IDs and is the exact input to cell relaxation.
    Each chemical element is therefore assigned independently with a global
    minimum-cost match under periodic boundary conditions.
    """
    by_element_start = {}
    for atom_id, element in start_elements.items():
        by_element_start.setdefault(element, []).append(atom_id)
    by_element_final = {}
    for index, element in enumerate(final_elements):
        by_element_final.setdefault(element, []).append(index)
    if set(by_element_start) != set(by_element_final):
        raise ValueError("element sets differ between atom and final cell dumps")

    cell_array = np.asarray(final_cell)
    coords_by_id = {}
    maximum_displacement = 0.0
    for element, atom_ids in by_element_start.items():
        final_indices = by_element_final[element]
        if len(atom_ids) != len(final_indices):
            raise ValueError(f"element count mismatch for {element}")
        initial = np.asarray([start_coords[atom_id] for atom_id in atom_ids])
        final = np.asarray([final_coords[index] for index in final_indices])
        delta = final[None, :, :] - initial[:, None, :]
        delta -= np.rint(delta)
        displacement = np.einsum("...j,jk->...k", delta, cell_array)
        costs = np.linalg.norm(displacement, axis=2)
        rows, columns = linear_sum_assignment(costs)
        assigned = costs[rows, columns]
        maximum_displacement = max(maximum_displacement, float(assigned.max(initial=0.0)))
        for row, column in zip(rows, columns):
            coords_by_id[atom_ids[row]] = tuple(final[column])
    return coords_by_id, maximum_displacement


def replace_data_coordinates(lines, atoms, coords_by_id, cell):
    if len(atoms) != len(coords_by_id):
        raise ValueError(f"Atom count mismatch: data={len(atoms)} dump={len(coords_by_id)}")
    output = list(lines)
    lx, ly, lz = cell[0][0], cell[1][1], cell[2][2]
    xy, xz, yz = cell[1][0], cell[2][0], cell[2][1]
    for i, line in enumerate(output):
        if "xlo xhi" in line:
            output[i] = f"{0.0:16.8f} {lx:16.8f} xlo xhi"
        elif "ylo yhi" in line:
            output[i] = f"{0.0:16.8f} {ly:16.8f} ylo yhi"
        elif "zlo zhi" in line:
            output[i] = f"{0.0:16.8f} {lz:16.8f} zlo zhi"
        elif "xy xz yz" in line:
            output[i] = f"{xy:16.8f} {xz:16.8f} {yz:16.8f} xy xz yz"
    for atom_id, (line_index, _) in atoms.items():
        x, y, z = cartesian(coords_by_id[atom_id], cell)
        values = output[line_index].split()
        values[4:7] = [f"{x:.8f}", f"{y:.8f}", f"{z:.8f}"]
        output[line_index] = " ".join(values)
    return "\n".join(output) + "\n"


def write_bonded_cif(path, name, atoms, atom_types, bonds, coords_by_id, cell):
    if len(atoms) != len(coords_by_id):
        raise ValueError(f"Atom count mismatch while writing {path}")
    a, b, c, alpha, beta, gamma = cell_parameters(cell)
    cif = [
        f"data_{name}_UFF4MOF_cell_bonded",
        "_audit_creation_method 'LAMMPS topology restored from data file and final cell dump'",
        "_symmetry_space_group_name_H-M 'P 1'",
        "_space_group_IT_number 1",
        f"_cell_length_a {a:.8f}",
        f"_cell_length_b {b:.8f}",
        f"_cell_length_c {c:.8f}",
        f"_cell_angle_alpha {alpha:.8f}",
        f"_cell_angle_beta {beta:.8f}",
        f"_cell_angle_gamma {gamma:.8f}",
        "loop_",
        "_space_group_symop_id",
        "_space_group_symop_operation_xyz",
        "1 x,y,z",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
    ]
    for atom_id in range(1, len(atoms) + 1):
        element = atom_types[atoms[atom_id][1]]
        frac = coords_by_id[atom_id]
        wrapped = tuple(value % 1.0 for value in frac)
        cif.append(f"{element}{atom_id} {element} {wrapped[0]:.8f} {wrapped[1]:.8f} {wrapped[2]:.8f}")
    cif += [
        "loop_",
        "_geom_bond_atom_site_label_1",
        "_geom_bond_atom_site_label_2",
        "_geom_bond_distance",
        "_geom_bond_site_symmetry_2",
        "_lammps_bond_type",
    ]
    for bond_type, atom1, atom2 in bonds:
        distance, translation = nearest_image_delta(coords_by_id[atom1], coords_by_id[atom2], cell)
        element1 = atom_types[atoms[atom1][1]]
        element2 = atom_types[atoms[atom2][1]]
        cif.append(
            f"{element1}{atom1} {element2}{atom2} "
            f"{distance:.6f} {periodic_code(translation)} {bond_type}"
        )
    path.write_text("\n".join(cif) + "\n")


def process_structure(workdir: Path, output_dir: Path):
    name = workdir.name
    data_path = workdir / f"data.{name}"
    dump_path = workdir / f"relaxed_{name}_cell_from_atom.lammpstrj"
    atom_dump_path = workdir / f"atom_relaxed_{name}_for_cell.lammpstrj"
    lines, atoms, atom_types, bonds = parse_data(data_path)
    elements, final_coords, cell, _ = final_frame(dump_path)
    start_coords, start_elements = parse_atom_relaxed_dump(atom_dump_path, atom_types)
    final_by_id, max_displacement = match_final_coordinates(
        start_coords, start_elements, elements, final_coords, cell
    )
    if max_displacement > 3.0:
        raise ValueError(f"unsafe coordinate assignment: maximum displacement {max_displacement:.3f} A")
    restored_data = workdir / f"{name}_UFF4MOF_cell.data"
    restored_data.write_text(replace_data_coordinates(lines, atoms, final_by_id, cell))
    write_bonded_cif(
        output_dir / f"{name}_UFF4MOF_cell_bonded.cif",
        name,
        atoms,
        atom_types,
        bonds,
        final_by_id,
        cell,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--metals", nargs="+", choices=METALS, default=METALS)
    args = parser.parse_args()
    total = succeeded = 0
    for metal in args.metals:
        root = ROOT / f"{metal}_3_lammps"
        output_dir = root / "All_cell_cifs"
        output_dir.mkdir(exist_ok=True)
        metal_total = metal_succeeded = 0
        for workdir in sorted(
            path for path in root.iterdir() if path.is_dir() and (path / f"data.{path.name}").exists()
        ):
            total += 1
            metal_total += 1
            try:
                process_structure(workdir, output_dir)
                succeeded += 1
                metal_succeeded += 1
            except Exception as exc:
                print(f"FAILED {metal} {workdir.name}: {exc}")
        print(f"{metal}: {metal_succeeded}/{metal_total} completed", flush=True)
    print(f"completed {succeeded}/{total}")


if __name__ == "__main__":
    main()
