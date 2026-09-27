"""Sparse empirical reference model and conformal-style rarity scoring."""

from __future__ import annotations

from dataclasses import dataclass, field
import gzip
import json
import math
from pathlib import Path
from statistics import median
from typing import Iterable

import numpy as np

from .geometry import Observation


@dataclass(frozen=True)
class HistogramSpec:
    minimum: float
    maximum: float
    width: float

    @property
    def n_bins(self) -> int:
        return int(math.ceil((self.maximum - self.minimum) / self.width))

    def bin_index(self, value: float) -> int | None:
        if not math.isfinite(value) or value < self.minimum or value > self.maximum:
            return None
        if value == self.maximum:
            return self.n_bins - 1
        index = int((value - self.minimum) // self.width)
        if 0 <= index < self.n_bins:
            return index
        return None

    def to_dict(self) -> dict:
        return {
            "minimum": self.minimum,
            "maximum": self.maximum,
            "width": self.width,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "HistogramSpec":
        return cls(float(data["minimum"]), float(data["maximum"]), float(data["width"]))


DEFAULT_SPECS = {
    "bond": HistogramSpec(0.5, 4.0, 0.02),
    "angle": HistogramSpec(0.0, 180.0, 1.0),
    "dihedral": HistogramSpec(0.0, 180.0, 2.0),
}


@dataclass
class SparseHistogram:
    total: int = 0
    counts: dict[int, int] = field(default_factory=dict)

    def add(self, index: int) -> None:
        self.total += 1
        self.counts[index] = self.counts.get(index, 0) + 1

    def rarity_p(self, index: int | None) -> float:
        """Finite-sample rarity probability for an observed histogram bin.

        The +1 correction prevents infinite anomaly values. An empty or
        out-of-range bin receives the finite-sample floor 1/(N+1).
        """

        if self.total <= 0:
            return 1.0
        if index is None:
            return 1.0 / (self.total + 1.0)
        observed = self.counts.get(index, 0)
        if observed == 0:
            return 1.0 / (self.total + 1.0)
        less_or_equal_mass = sum(count for count in self.counts.values() if count <= observed)
        return min(1.0, (1.0 + less_or_equal_mass) / (self.total + 1.0))

    def to_dict(self) -> dict:
        return {
            "total": self.total,
            "counts": {str(k): v for k, v in sorted(self.counts.items())},
        }

    @classmethod
    def from_dict(cls, data: dict) -> "SparseHistogram":
        return cls(
            total=int(data["total"]),
            counts={int(k): int(v) for k, v in data["counts"].items()},
        )


def _quantile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    return float(np.quantile(np.asarray(values, dtype=float), q, method="linear"))


@dataclass
class ReferenceModel:
    """Hierarchical reference histograms for MOF internal coordinates."""

    specs: dict[str, HistogramSpec] = field(
        default_factory=lambda: dict(DEFAULT_SPECS)
    )
    min_count: int = 100
    min_coverage: float = 0.50
    min_kind_coverage: float = 0.25
    max_environment_level: int = 1
    structure_quantile: float = 0.99
    bond_scale: float = 1.15
    prune_below: int = 5
    tables: dict[str, list[dict[str, SparseHistogram]]] = field(default_factory=dict)
    threshold: float | None = None
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.tables:
            self.tables = {kind: [dict(), dict(), dict()] for kind in self.specs}
        if self.min_count < 1:
            raise ValueError("min_count must be at least 1")
        if not (0.0 <= self.min_coverage <= 1.0):
            raise ValueError("min_coverage must lie in [0, 1]")
        if not (0.0 <= self.min_kind_coverage <= 1.0):
            raise ValueError("min_kind_coverage must lie in [0, 1]")
        if self.max_environment_level not in (0, 1, 2):
            raise ValueError("max_environment_level must be 0, 1, or 2")
        if not (0.0 < self.structure_quantile <= 1.0):
            raise ValueError("structure_quantile must lie in (0, 1]")

    def add_observations(self, observations: Iterable[Observation]) -> None:
        for observation in observations:
            if observation.kind not in self.specs:
                continue
            index = self.specs[observation.kind].bin_index(observation.value)
            if index is None:
                continue
            levels = self.tables[observation.kind]
            for level, key in enumerate(observation.keys[: len(levels)]):
                histogram = levels[level].setdefault(key, SparseHistogram())
                histogram.add(index)

    def prune(self) -> dict[str, int]:
        """Remove environment histograms too small to be useful after fitting."""

        removed: dict[str, int] = {}
        for kind, levels in self.tables.items():
            count = 0
            for table in levels:
                doomed = [key for key, hist in table.items() if hist.total < self.prune_below]
                for key in doomed:
                    del table[key]
                count += len(doomed)
            removed[kind] = count
        return removed

    def _score_one(self, observation: Observation) -> dict:
        spec = self.specs[observation.kind]
        index = spec.bin_index(observation.value)
        selected = None
        for level, key in enumerate(observation.keys):
            if level > self.max_environment_level:
                break
            if level >= len(self.tables[observation.kind]):
                break
            histogram = self.tables[observation.kind][level].get(key)
            if histogram is not None and histogram.total >= self.min_count:
                selected = (level, key, histogram)
                break

        base = {
            "kind": observation.kind,
            "value": float(observation.value),
            "atoms": list(observation.atoms),
        }
        if selected is None:
            return {
                **base,
                "supported": False,
                "level": None,
                "key": None,
                "reference_count": 0,
                "p_value": None,
                "anomaly": None,
            }

        level, key, histogram = selected
        p_value = histogram.rarity_p(index)
        anomaly = -math.log10(max(p_value, np.finfo(float).tiny))
        return {
            **base,
            "supported": True,
            "level": int(level),
            "key": key,
            "reference_count": int(histogram.total),
            "p_value": float(p_value),
            "anomaly": float(anomaly),
        }

    def score_observations(
        self,
        observations: Iterable[Observation],
        *,
        structure_id: str | None = None,
        top_n: int = 20,
    ) -> dict:
        local = [self._score_one(observation) for observation in observations]
        supported = [item for item in local if item["supported"]]
        total = len(local)
        coverage = len(supported) / total if total else 0.0

        kind_summary: dict[str, dict] = {}
        kind_scores: list[float] = []
        for kind in self.specs:
            all_kind = [item for item in local if item["kind"] == kind]
            supported_kind = [item for item in all_kind if item["supported"]]
            anomalies = [float(item["anomaly"]) for item in supported_kind]
            q_score = _quantile(anomalies, self.structure_quantile)
            if q_score is not None:
                kind_scores.append(q_score)
            kind_summary[kind] = {
                "n_total": len(all_kind),
                "n_supported": len(supported_kind),
                "coverage": len(supported_kind) / len(all_kind) if all_kind else 0.0,
                "median": float(median(anomalies)) if anomalies else None,
                "q95": _quantile(anomalies, 0.95),
                "q99": _quantile(anomalies, 0.99),
                "maximum": max(anomalies) if anomalies else None,
                "structure_quantile": q_score,
            }

        score = float(np.mean(kind_scores)) if kind_scores else None
        populated_kind_coverages = [
            summary["coverage"]
            for summary in kind_summary.values()
            if summary["n_total"] > 0
        ]
        kind_domain_failure = any(
            value < self.min_kind_coverage for value in populated_kind_coverages
        )
        if coverage < self.min_coverage or kind_domain_failure or score is None:
            status = "out_of_domain"
        elif self.threshold is None:
            status = "unclassified"
        elif score > self.threshold:
            status = "geometric_outlier"
        else:
            status = "in_distribution"

        ranked = sorted(
            supported,
            key=lambda item: float(item["anomaly"]),
            reverse=True,
        )
        return {
            "schema_version": 1,
            "structure_id": structure_id,
            "status": status,
            "score": score,
            "threshold": self.threshold,
            "coverage": coverage,
            "n_observations": total,
            "n_supported": len(supported),
            "min_reference_count": self.min_count,
            "min_kind_coverage": self.min_kind_coverage,
            "max_environment_level": self.max_environment_level,
            "structure_quantile": self.structure_quantile,
            "kind_summary": kind_summary,
            "top_anomalies": ranked[: max(0, top_n)],
            "unsupported_examples": [item for item in local if not item["supported"]][
                : max(0, min(top_n, 10))
            ],
            "model_metadata": self.metadata,
        }

    def calibrate(self, structure_scores: Iterable[float], quantile: float = 0.99) -> float:
        scores = [float(value) for value in structure_scores if math.isfinite(float(value))]
        if len(scores) < 5:
            raise ValueError("at least five finite calibration scores are required")
        if not (0.0 < quantile < 1.0):
            raise ValueError("calibration quantile must lie in (0, 1)")
        self.threshold = float(np.quantile(scores, quantile, method="higher"))
        self.metadata["calibration"] = {
            "n_structures": len(scores),
            "quantile": quantile,
            "threshold": self.threshold,
        }
        return self.threshold

    def to_dict(self) -> dict:
        return {
            "format": "mofrelax-reference-model",
            "format_version": 1,
            "specs": {kind: spec.to_dict() for kind, spec in self.specs.items()},
            "min_count": self.min_count,
            "min_coverage": self.min_coverage,
            "min_kind_coverage": self.min_kind_coverage,
            "max_environment_level": self.max_environment_level,
            "structure_quantile": self.structure_quantile,
            "bond_scale": self.bond_scale,
            "prune_below": self.prune_below,
            "threshold": self.threshold,
            "metadata": self.metadata,
            "tables": {
                kind: [
                    {key: histogram.to_dict() for key, histogram in sorted(table.items())}
                    for table in levels
                ]
                for kind, levels in self.tables.items()
            },
        }

    @classmethod
    def from_dict(cls, data: dict) -> "ReferenceModel":
        if data.get("format") != "mofrelax-reference-model":
            raise ValueError("not a mofrelax reference model")
        tables = {
            kind: [
                {key: SparseHistogram.from_dict(hist) for key, hist in table.items()}
                for table in levels
            ]
            for kind, levels in data["tables"].items()
        }
        return cls(
            specs={kind: HistogramSpec.from_dict(spec) for kind, spec in data["specs"].items()},
            min_count=int(data["min_count"]),
            min_coverage=float(data["min_coverage"]),
            min_kind_coverage=float(data.get("min_kind_coverage", 0.25)),
            max_environment_level=int(data.get("max_environment_level", 1)),
            structure_quantile=float(data["structure_quantile"]),
            bond_scale=float(data["bond_scale"]),
            prune_below=int(data.get("prune_below", 5)),
            tables=tables,
            threshold=float(data["threshold"]) if data.get("threshold") is not None else None,
            metadata=dict(data.get("metadata", {})),
        )

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(temporary, "wt", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, separators=(",", ":"), sort_keys=True)
        temporary.replace(path)

    @classmethod
    def load(cls, path: str | Path) -> "ReferenceModel":
        path = Path(path)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as handle:
            return cls.from_dict(json.load(handle))
