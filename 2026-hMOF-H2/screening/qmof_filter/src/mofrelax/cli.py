"""Command-line interface for QMOF data preparation, training, and scoring."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

from .io import ManifestRow, find_structure_files, load_manifest, make_manifest
from .model import ReferenceModel
from .pipeline import score_rows, train_reference_model, write_diagnostics
from .qmof import download_files, list_files


def _print_json(data) -> None:
    print(json.dumps(data, indent=2, sort_keys=True))


def cmd_qmof_files(args) -> int:
    files = list_files()
    _print_json([{k: v for k, v in item.items() if k != "download_url"} for item in files])
    return 0


def cmd_qmof_download(args) -> int:
    result = download_files(
        args.output_dir,
        names=set(args.name) if args.name else None,
        match=args.match,
        extract=args.extract,
    )
    _print_json(result)
    return 0


def cmd_make_manifest(args) -> int:
    result = make_manifest(args.cif_dir, args.output, metadata=args.metadata)
    _print_json(result)
    return 0


def cmd_train(args) -> int:
    sources = set(args.sources) if args.sources else None
    rows = load_manifest(args.manifest, sources=sources)
    if not rows:
        raise ValueError("no manifest rows remained after include/source filtering")
    metadata = {
        "reference_dataset": "QMOF",
        "qmof_version": args.qmof_version,
        "source_filter": sorted(sources) if sources else None,
        "manifest": str(Path(args.manifest).resolve()),
    }
    model, diagnostics = train_reference_model(
        rows,
        min_count=args.min_count,
        min_coverage=args.min_coverage,
        min_kind_coverage=args.min_kind_coverage,
        max_environment_level=args.max_environment_level,
        structure_quantile=args.structure_quantile,
        calibration_fraction=args.calibration_fraction,
        calibration_quantile=args.calibration_quantile,
        bond_scale=args.bond_scale,
        prune_below=args.prune_below,
        workers=args.workers,
        metadata=metadata,
    )
    model.save(args.model)
    diagnostics_path = args.diagnostics or str(Path(args.model).with_suffix(".diagnostics.json"))
    write_diagnostics(diagnostics_path, diagnostics)
    _print_json(
        {
            "model": str(args.model),
            "diagnostics": diagnostics_path,
            "threshold": model.threshold,
            "n_train_successful": diagnostics["n_train_successful"],
            "n_calibration_supported": diagnostics["n_calibration_supported"],
            "n_observations": diagnostics["n_observations"],
            "warning": diagnostics["warning"],
        }
    )
    return 0


def _rows_from_input(path: str | Path) -> list[ManifestRow]:
    files = find_structure_files(path)
    return [ManifestRow(item.stem, item.resolve(), "", item.stem, True) for item in files]


def _write_reports(path: str | Path, reports: list[dict]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if len(reports) == 1 and path.suffix.lower() != ".jsonl":
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(reports[0], handle, indent=2, sort_keys=True)
    else:
        with open(path, "w", encoding="utf-8") as handle:
            for report in reports:
                handle.write(json.dumps(report, sort_keys=True) + "\n")


def _write_summary(path: str | Path, reports: list[dict]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "structure_id",
        "path",
        "status",
        "score",
        "threshold",
        "coverage",
        "n_observations",
        "n_supported",
        "error",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for report in reports:
            writer.writerow({key: report.get(key) for key in fields})


def cmd_score(args) -> int:
    model = ReferenceModel.load(args.model)
    rows = _rows_from_input(args.input)
    if not rows:
        raise ValueError("no CIF files found")
    reports = score_rows(rows, model, workers=args.workers, top_n=args.top)
    _write_reports(args.output, reports)
    if args.summary_csv:
        _write_summary(args.summary_csv, reports)
    counts: dict[str, int] = {}
    for report in reports:
        counts[report["status"]] = counts.get(report["status"], 0) + 1
    _print_json(
        {
            "output": str(args.output),
            "summary_csv": args.summary_csv,
            "n_structures": len(reports),
            "status_counts": counts,
        }
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mofrelax",
        description="QMOF-calibrated geometric plausibility filters for MOFs",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    sub = subparsers.add_parser("qmof-files", help="list files in the official QMOF Figshare article")
    sub.set_defaults(func=cmd_qmof_files)

    sub = subparsers.add_parser("qmof-download", help="download selected official QMOF files")
    sub.add_argument("output_dir")
    sub.add_argument("--name", action="append", default=[], help="exact Figshare filename; repeatable")
    sub.add_argument("--match", help="regular expression matched against Figshare filenames")
    sub.add_argument("--extract", action="store_true", help="safely extract downloaded archives")
    sub.set_defaults(func=cmd_qmof_download)

    sub = subparsers.add_parser("make-manifest", help="index QMOF CIFs and optional JSON metadata")
    sub.add_argument("--cif-dir", required=True)
    sub.add_argument("--metadata")
    sub.add_argument("--output", required=True)
    sub.set_defaults(func=cmd_make_manifest)

    sub = subparsers.add_parser("train", help="fit and calibrate a geometry reference model")
    sub.add_argument("--manifest", required=True)
    sub.add_argument("--model", required=True)
    sub.add_argument("--diagnostics")
    sub.add_argument("--sources", nargs="+", default=[])
    sub.add_argument("--qmof-version", default="unspecified")
    sub.add_argument("--min-count", type=int, default=100)
    sub.add_argument("--min-coverage", type=float, default=0.50)
    sub.add_argument("--min-kind-coverage", type=float, default=0.25)
    sub.add_argument("--max-environment-level", type=int, choices=(0, 1, 2), default=1)
    sub.add_argument("--structure-quantile", type=float, default=0.99)
    sub.add_argument("--calibration-fraction", type=float, default=0.10)
    sub.add_argument("--calibration-quantile", type=float, default=0.99)
    sub.add_argument("--bond-scale", type=float, default=1.15)
    sub.add_argument("--prune-below", type=int, default=5)
    sub.add_argument("--workers", type=int, default=1)
    sub.set_defaults(func=cmd_train)

    sub = subparsers.add_parser("score", help="score one CIF or a directory of CIFs")
    sub.add_argument("--model", required=True)
    sub.add_argument("--input", required=True)
    sub.add_argument("--output", required=True)
    sub.add_argument("--summary-csv")
    sub.add_argument("--top", type=int, default=20)
    sub.add_argument("--workers", type=int, default=1)
    sub.set_defaults(func=cmd_score)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except KeyboardInterrupt:
        return 130
    except Exception as exc:
        print(f"mofrelax: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
