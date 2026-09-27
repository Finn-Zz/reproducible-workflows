"""Training and scoring workflows, including process-level parallelism."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
from typing import Iterable, Iterator

from .geometry import Observation
from .io import ManifestRow, observations_from_cif
from .model import ReferenceModel


@dataclass(frozen=True)
class ExtractionResult:
    structure_id: str
    path: str
    observations: tuple[Observation, ...]
    error: str | None = None


def _extract_worker(task: tuple[str, str, float]) -> ExtractionResult:
    structure_id, path, bond_scale = task
    try:
        observations = tuple(observations_from_cif(path, bond_scale=bond_scale))
        return ExtractionResult(structure_id, path, observations, None)
    except Exception as exc:  # batch workflows must record and continue
        return ExtractionResult(structure_id, path, tuple(), f"{type(exc).__name__}: {exc}")


def extract_many(
    rows: Iterable[ManifestRow],
    *,
    bond_scale: float,
    workers: int = 1,
) -> Iterator[ExtractionResult]:
    tasks = [(row.structure_id, str(row.path), bond_scale) for row in rows]
    if workers <= 1:
        for task in tasks:
            yield _extract_worker(task)
        return
    with ProcessPoolExecutor(max_workers=workers) as executor:
        yield from executor.map(_extract_worker, tasks, chunksize=1)


def stable_fraction(value: str) -> float:
    digest = hashlib.sha256(value.encode("utf-8")).digest()
    integer = int.from_bytes(digest[:8], "big")
    return integer / float(2**64)


def split_rows(
    rows: Iterable[ManifestRow], calibration_fraction: float
) -> tuple[list[ManifestRow], list[ManifestRow]]:
    if not (0.0 < calibration_fraction < 1.0):
        raise ValueError("calibration_fraction must lie in (0, 1)")
    train: list[ManifestRow] = []
    calibration: list[ManifestRow] = []
    for row in rows:
        group = row.split_group or row.structure_id
        if stable_fraction(group) < calibration_fraction:
            calibration.append(row)
        else:
            train.append(row)
    if not train or not calibration:
        raise ValueError("deterministic split produced an empty training or calibration set")
    return train, calibration


def train_reference_model(
    rows: list[ManifestRow],
    *,
    min_count: int = 100,
    min_coverage: float = 0.50,
    min_kind_coverage: float = 0.25,
    max_environment_level: int = 1,
    structure_quantile: float = 0.99,
    calibration_fraction: float = 0.10,
    calibration_quantile: float = 0.99,
    bond_scale: float = 1.15,
    prune_below: int = 5,
    workers: int = 1,
    metadata: dict | None = None,
) -> tuple[ReferenceModel, dict]:
    train_rows, calibration_rows = split_rows(rows, calibration_fraction)
    model = ReferenceModel(
        min_count=min_count,
        min_coverage=min_coverage,
        min_kind_coverage=min_kind_coverage,
        max_environment_level=max_environment_level,
        structure_quantile=structure_quantile,
        bond_scale=bond_scale,
        prune_below=prune_below,
        metadata=dict(metadata or {}),
    )
    model.metadata.update(
        {
            "created_at": datetime.now(timezone.utc).isoformat(),
            "software": "qmof-relax-filter 0.1.0",
            "python": platform.python_version(),
            "reference_structure_count_requested": len(train_rows),
            "calibration_structure_count_requested": len(calibration_rows),
            "bond_scale": bond_scale,
            "max_environment_level": max_environment_level,
            "scientific_scope": (
                "Geometric reference rarity only; not an energy, synthesizability, "
                "thermodynamic-stability, or dynamic-stability prediction."
            ),
        }
    )

    train_errors: list[dict] = []
    n_train_success = 0
    n_observations = 0
    for result in extract_many(train_rows, bond_scale=bond_scale, workers=workers):
        if result.error:
            train_errors.append({"structure_id": result.structure_id, "path": result.path, "error": result.error})
            continue
        model.add_observations(result.observations)
        n_train_success += 1
        n_observations += len(result.observations)
    removed = model.prune()

    calibration_errors: list[dict] = []
    calibration_scores: list[float] = []
    calibration_reports: list[dict] = []
    for result in extract_many(calibration_rows, bond_scale=bond_scale, workers=workers):
        if result.error:
            calibration_errors.append(
                {"structure_id": result.structure_id, "path": result.path, "error": result.error}
            )
            continue
        report = model.score_observations(
            result.observations, structure_id=result.structure_id, top_n=0
        )
        if report["score"] is not None and report["coverage"] >= model.min_coverage:
            calibration_scores.append(float(report["score"]))
            calibration_reports.append(
                {
                    "structure_id": result.structure_id,
                    "score": report["score"],
                    "coverage": report["coverage"],
                }
            )

    warning = None
    if len(calibration_scores) >= 5:
        model.calibrate(calibration_scores, quantile=calibration_quantile)
    else:
        warning = (
            "Fewer than five supported calibration structures were available; "
            "the model was saved without a structure-level threshold."
        )
        model.metadata["calibration"] = {
            "n_structures": len(calibration_scores),
            "quantile": calibration_quantile,
            "threshold": None,
            "warning": warning,
        }

    model.metadata.update(
        {
            "reference_structure_count_successful": n_train_success,
            "reference_observation_count": n_observations,
            "reference_extraction_failures": len(train_errors),
            "calibration_extraction_failures": len(calibration_errors),
            "pruned_environment_count": removed,
        }
    )
    diagnostics = {
        "n_train_requested": len(train_rows),
        "n_train_successful": n_train_success,
        "n_calibration_requested": len(calibration_rows),
        "n_calibration_supported": len(calibration_scores),
        "n_observations": n_observations,
        "threshold": model.threshold,
        "warning": warning,
        "train_errors": train_errors,
        "calibration_errors": calibration_errors,
        "calibration_reports": calibration_reports,
        "pruned_environment_count": removed,
    }
    return model, diagnostics


def score_cif(
    path: str | Path,
    model: ReferenceModel,
    *,
    top_n: int = 20,
) -> dict:
    path = Path(path)
    observations = observations_from_cif(path, bond_scale=model.bond_scale)
    report = model.score_observations(observations, structure_id=path.stem, top_n=top_n)
    report["path"] = str(path)
    return report


def score_rows(
    rows: Iterable[ManifestRow],
    model: ReferenceModel,
    *,
    workers: int = 1,
    top_n: int = 20,
) -> list[dict]:
    reports: list[dict] = []
    for result in extract_many(rows, bond_scale=model.bond_scale, workers=workers):
        if result.error:
            reports.append(
                {
                    "schema_version": 1,
                    "structure_id": result.structure_id,
                    "path": result.path,
                    "status": "error",
                    "error": result.error,
                    "score": None,
                    "coverage": 0.0,
                }
            )
            continue
        report = model.score_observations(
            result.observations, structure_id=result.structure_id, top_n=top_n
        )
        report["path"] = result.path
        reports.append(report)
    return reports


def write_diagnostics(path: str | Path, diagnostics: dict) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(diagnostics, handle, indent=2, sort_keys=True)
