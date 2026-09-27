"""Periodic graph and internal-coordinate extraction.

The numerical geometry functions are library-independent. ASE is imported only
by :func:`infer_edges_ase`, which keeps the model and tests lightweight.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import math
from typing import Iterable, Sequence

import numpy as np

from .elements import element_class, is_metal


@dataclass(frozen=True)
class StructureData:
    """Minimal periodic structure representation.

    ``cell`` uses ASE's convention: lattice vectors are rows of a 3x3 matrix.
    """

    symbols: tuple[str, ...]
    positions: np.ndarray
    cell: np.ndarray
    pbc: tuple[bool, bool, bool] = (True, True, True)
    identifier: str | None = None

    def __post_init__(self) -> None:
        positions = np.asarray(self.positions, dtype=float)
        cell = np.asarray(self.cell, dtype=float)
        if positions.shape != (len(self.symbols), 3):
            raise ValueError("positions must have shape (n_atoms, 3)")
        if cell.shape != (3, 3):
            raise ValueError("cell must have shape (3, 3)")
        if abs(float(np.linalg.det(cell))) < 1.0e-12:
            raise ValueError("cell is singular")
        object.__setattr__(self, "positions", positions)
        object.__setattr__(self, "cell", cell)


@dataclass(frozen=True)
class Observation:
    kind: str
    value: float
    atoms: tuple[int, ...]
    keys: tuple[str, ...]  # fine -> coarse

    def to_dict(self) -> dict:
        return {
            "kind": self.kind,
            "value": float(self.value),
            "atoms": list(self.atoms),
            "keys": list(self.keys),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Observation":
        return cls(
            kind=str(data["kind"]),
            value=float(data["value"]),
            atoms=tuple(int(x) for x in data["atoms"]),
            keys=tuple(str(x) for x in data["keys"]),
        )


def minimum_image_vector(
    delta: np.ndarray, cell: np.ndarray, pbc: Sequence[bool]
) -> np.ndarray:
    """Return the minimum-image displacement for a triclinic periodic cell."""

    delta = np.asarray(delta, dtype=float)
    fractional = delta @ np.linalg.inv(cell)
    for axis, periodic in enumerate(pbc):
        if periodic:
            fractional[axis] -= np.rint(fractional[axis])
    return fractional @ cell


def bond_vector(structure: StructureData, i: int, j: int) -> np.ndarray:
    return minimum_image_vector(
        structure.positions[j] - structure.positions[i],
        structure.cell,
        structure.pbc,
    )


def distance(structure: StructureData, i: int, j: int) -> float:
    return float(np.linalg.norm(bond_vector(structure, i, j)))


def angle_degrees(structure: StructureData, i: int, j: int, k: int) -> float:
    a = bond_vector(structure, j, i)
    b = bond_vector(structure, j, k)
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 1.0e-14:
        return float("nan")
    cosine = float(np.dot(a, b) / denom)
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def dihedral_degrees(
    structure: StructureData, i: int, j: int, k: int, l: int
) -> float:
    """Return an unsigned, folded torsion in [0, 180] degrees."""

    p1 = structure.positions[j]
    p0 = p1 + bond_vector(structure, j, i)
    p2 = p1 + bond_vector(structure, j, k)
    p3 = p2 + bond_vector(structure, k, l)

    b0 = p1 - p0
    b1 = p2 - p1
    b2 = p3 - p2
    norm_b1 = float(np.linalg.norm(b1))
    if norm_b1 <= 1.0e-14:
        return float("nan")
    unit_b1 = b1 / norm_b1
    v = b0 - np.dot(b0, unit_b1) * unit_b1
    w = b2 - np.dot(b2, unit_b1) * unit_b1
    nv = float(np.linalg.norm(v))
    nw = float(np.linalg.norm(w))
    if nv <= 1.0e-14 or nw <= 1.0e-14:
        return float("nan")
    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(unit_b1, v), w))
    signed = math.degrees(math.atan2(y, x))
    return abs(float(signed))


def adjacency_from_edges(n_atoms: int, edges: Iterable[tuple[int, int]]) -> list[set[int]]:
    adjacency = [set() for _ in range(n_atoms)]
    for i, j in edges:
        if i == j:
            continue
        if not (0 <= i < n_atoms and 0 <= j < n_atoms):
            raise IndexError(f"edge {(i, j)} is outside a {n_atoms}-atom structure")
        adjacency[i].add(j)
        adjacency[j].add(i)
    return adjacency


def infer_edges_ase(atoms, bond_scale: float = 1.15) -> list[tuple[int, int]]:
    """Infer periodic covalent edges with ASE's natural cutoffs.

    The scale is intentionally configurable and must be recorded with a model.
    A sensitivity analysis over plausible scales is recommended for publication.
    """

    try:
        from ase.neighborlist import natural_cutoffs, neighbor_list
    except ImportError as exc:  # pragma: no cover - exercised in integration
        raise ImportError("CIF processing requires ASE: pip install ase") from exc

    if bond_scale <= 0:
        raise ValueError("bond_scale must be positive")
    cutoffs = natural_cutoffs(atoms, mult=bond_scale)
    i_array, j_array = neighbor_list("ij", atoms, cutoffs)
    edges = {tuple(sorted((int(i), int(j)))) for i, j in zip(i_array, j_array) if i != j}
    return sorted(edges)


def _atom_labels(symbol: str, coordination: int) -> tuple[str, str, str]:
    return (
        f"{symbol}|cn{coordination}",
        symbol,
        element_class(symbol),
    )


def _canonical(sequence: Sequence[str]) -> str:
    forward = "-".join(sequence)
    reverse = "-".join(reversed(sequence))
    return min(forward, reverse)


def _bond_keys(symbols: Sequence[str], degrees: Sequence[int], i: int, j: int) -> tuple[str, ...]:
    li = _atom_labels(symbols[i], degrees[i])
    lj = _atom_labels(symbols[j], degrees[j])
    return tuple(_canonical((li[level], lj[level])) for level in range(3))


def _angle_keys(
    symbols: Sequence[str], degrees: Sequence[int], i: int, j: int, k: int
) -> tuple[str, ...]:
    labels = [_atom_labels(symbols[x], degrees[x]) for x in (i, j, k)]
    keys = []
    for level in range(3):
        outer = sorted((labels[0][level], labels[2][level]))
        keys.append(f"{outer[0]}-{labels[1][level]}-{outer[1]}")
    return tuple(keys)


def _dihedral_keys(
    symbols: Sequence[str], degrees: Sequence[int], path: tuple[int, int, int, int]
) -> tuple[str, ...]:
    labels = [_atom_labels(symbols[x], degrees[x]) for x in path]
    return tuple(_canonical([label[level] for label in labels]) for level in range(3))


def extract_observations(
    structure: StructureData,
    edges: Iterable[tuple[int, int]],
    *,
    include_bonds: bool = True,
    include_angles: bool = True,
    include_dihedrals: bool = True,
    include_metal_centered_dihedrals: bool = False,
) -> list[Observation]:
    """Extract unique internal coordinates and hierarchical environment keys."""

    edge_set = {tuple(sorted(edge)) for edge in edges if edge[0] != edge[1]}
    adjacency = adjacency_from_edges(len(structure.symbols), edge_set)
    degrees = [len(neighbors) for neighbors in adjacency]
    symbols = structure.symbols
    observations: list[Observation] = []

    if include_bonds:
        for i, j in sorted(edge_set):
            value = distance(structure, i, j)
            if math.isfinite(value):
                observations.append(
                    Observation("bond", value, (i, j), _bond_keys(symbols, degrees, i, j))
                )

    if include_angles:
        for j, neighbors in enumerate(adjacency):
            for i, k in combinations(sorted(neighbors), 2):
                value = angle_degrees(structure, i, j, k)
                if math.isfinite(value):
                    observations.append(
                        Observation(
                            "angle", value, (i, j, k), _angle_keys(symbols, degrees, i, j, k)
                        )
                    )

    if include_dihedrals:
        seen_paths: set[tuple[int, int, int, int]] = set()
        for j, k in sorted(edge_set):
            if not include_metal_centered_dihedrals and (
                is_metal(symbols[j]) or is_metal(symbols[k])
            ):
                continue
            for i in adjacency[j] - {k}:
                for l in adjacency[k] - {j}:
                    if i == l:
                        continue
                    path = (i, j, k, l)
                    canonical = min(path, tuple(reversed(path)))
                    if canonical in seen_paths:
                        continue
                    seen_paths.add(canonical)
                    value = dihedral_degrees(structure, *path)
                    if math.isfinite(value):
                        observations.append(
                            Observation(
                                "dihedral",
                                value,
                                path,
                                _dihedral_keys(symbols, degrees, path),
                            )
                        )
    return observations
