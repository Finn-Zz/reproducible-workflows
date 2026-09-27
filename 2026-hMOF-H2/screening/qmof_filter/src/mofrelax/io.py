"""Structure I/O and QMOF manifest construction."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np

from .geometry import StructureData, extract_observations, infer_edges_ase


CIF_SUFFIXES = {".cif", ".mcif"}


def read_structure(path: str | Path) -> tuple[StructureData, object]:
    """Read a structure through ASE and return both neutral and ASE forms."""

    try:
        from ase.io import read
    except ImportError as exc:  # pragma: no cover - integration only
        raise ImportError("CIF processing requires ASE: pip install ase") from exc

    path = Path(path)
    atoms = read(str(path), index=0)
    if len(atoms) == 0:
        raise ValueError(f"no atoms found in {path}")
    cell = np.asarray(atoms.cell.array, dtype=float)
    if abs(float(np.linalg.det(cell))) < 1.0e-12:
        raise ValueError(f"structure has no valid 3D unit cell: {path}")
    data = StructureData(
        symbols=tuple(atoms.get_chemical_symbols()),
        positions=np.asarray(atoms.positions, dtype=float),
        cell=cell,
        pbc=tuple(bool(x) for x in atoms.pbc),
        identifier=path.stem,
    )
    return data, atoms


def observations_from_cif(path: str | Path, bond_scale: float = 1.15):
    structure, atoms = read_structure(path)
    edges = infer_edges_ase(atoms, bond_scale=bond_scale)
    return extract_observations(structure, edges)


def find_structure_files(path: str | Path, recursive: bool = True) -> list[Path]:
    path = Path(path)
    if path.is_file():
        return [path]
    if not path.is_dir():
        raise FileNotFoundError(path)
    iterator = path.rglob("*") if recursive else path.glob("*")
    return sorted(item for item in iterator if item.is_file() and item.suffix.lower() in CIF_SUFFIXES)


def _walk_scalars(data, prefix: str = "") -> Iterator[tuple[str, object]]:
    if isinstance(data, dict):
        for key, value in data.items():
            new_prefix = f"{prefix}.{key}" if prefix else str(key)
            yield from _walk_scalars(value, new_prefix)
    elif isinstance(data, (str, int, float, bool)) or data is None:
        yield prefix, data


def _lookup_scalar(record: dict, candidates: Iterable[str]):
    values: dict[str, object] = {}
    for key, value in _walk_scalars(record):
        leaf = key.rsplit(".", 1)[-1].lower()
        if value not in (None, "") and leaf not in values:
            values[leaf] = value
    for candidate in candidates:
        if candidate.lower() in values:
            return values[candidate.lower()]
    return None


def load_qmof_records(path: str | Path) -> list[dict]:
    """Load common QMOF JSON layouts without hard-coding one release schema."""

    with open(path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    records: list[dict] = []
    if isinstance(data, list):
        records = [dict(item) for item in data if isinstance(item, dict)]
    elif isinstance(data, dict):
        if isinstance(data.get("data"), list):
            records = [dict(item) for item in data["data"] if isinstance(item, dict)]
        else:
            for key, value in data.items():
                if isinstance(value, dict):
                    record = dict(value)
                    record.setdefault("_record_key", key)
                    records.append(record)
    else:
        raise ValueError("unsupported QMOF metadata JSON layout")
    if not records:
        raise ValueError(f"no metadata records found in {path}")
    return records


def _record_identifiers(record: dict) -> set[str]:
    values = {
        _lookup_scalar(record, ("qmof_id", "id")),
        _lookup_scalar(record, ("name", "refcode", "entry", "filename", "file_name")),
        record.get("_record_key"),
    }
    identifiers: set[str] = set()
    for value in values:
        if value in (None, ""):
            continue
        text = Path(str(value)).stem
        identifiers.add(text)
        identifiers.add(text.lower())
    return identifiers


@dataclass(frozen=True)
class ManifestRow:
    structure_id: str
    path: Path
    source: str = ""
    split_group: str = ""
    include: bool = True


def make_manifest(
    cif_dir: str | Path,
    output: str | Path,
    metadata: str | Path | None = None,
) -> dict:
    files = find_structure_files(cif_dir)
    lookup: dict[str, list[dict]] = {}
    if metadata is not None:
        for record in load_qmof_records(metadata):
            for identifier in _record_identifiers(record):
                lookup.setdefault(identifier, []).append(record)

    rows: list[ManifestRow] = []
    matched = ambiguous = 0
    for path in files:
        candidates = lookup.get(path.stem, []) or lookup.get(path.stem.lower(), [])
        record = None
        unique_candidates = {id(item): item for item in candidates}
        if len(unique_candidates) == 1:
            record = next(iter(unique_candidates.values()))
            matched += 1
        elif len(unique_candidates) > 1:
            ambiguous += 1

        if record is None:
            structure_id = path.stem
            source = ""
        else:
            structure_id = str(
                _lookup_scalar(record, ("qmof_id", "id", "name", "refcode")) or path.stem
            )
            source = str(_lookup_scalar(record, ("source", "data_source", "database")) or "")
        rows.append(
            ManifestRow(
                structure_id=structure_id,
                path=path.resolve(),
                source=source,
                split_group=structure_id,
                include=True,
            )
        )

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["structure_id", "path", "source", "split_group", "include"],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "structure_id": row.structure_id,
                    "path": str(row.path),
                    "source": row.source,
                    "split_group": row.split_group,
                    "include": "1" if row.include else "0",
                }
            )
    return {
        "n_files": len(files),
        "n_metadata_matched": matched,
        "n_metadata_ambiguous": ambiguous,
        "output": str(output),
    }


def _truthy(value: str | None) -> bool:
    return str(value or "").strip().lower() not in {"", "0", "false", "no", "exclude"}


def load_manifest(
    path: str | Path,
    *,
    sources: set[str] | None = None,
) -> list[ManifestRow]:
    path = Path(path)
    rows: list[ManifestRow] = []
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"structure_id", "path"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise ValueError(f"manifest must contain columns {sorted(required)}")
        for raw in reader:
            include = _truthy(raw.get("include", "1"))
            source = str(raw.get("source", "")).strip()
            if not include or (sources is not None and source not in sources):
                continue
            item_path = Path(str(raw["path"]))
            if not item_path.is_absolute():
                item_path = (path.parent / item_path).resolve()
            structure_id = str(raw["structure_id"]).strip()
            split_group = str(raw.get("split_group") or structure_id).strip()
            rows.append(
                ManifestRow(structure_id, item_path, source, split_group, include=True)
            )
    return rows
