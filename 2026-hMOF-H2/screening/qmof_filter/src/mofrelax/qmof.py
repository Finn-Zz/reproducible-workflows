"""Minimal Figshare client for the separately distributed QMOF dataset."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import tarfile
import urllib.request
import zipfile


FIGSHARE_ARTICLE_ID = 13147324
FIGSHARE_API = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE_ID}"


def list_files() -> list[dict]:
    request = urllib.request.Request(
        FIGSHARE_API,
        headers={"User-Agent": "qmof-relax-filter/0.1 (+academic software)"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = json.load(response)
    files = []
    for item in payload.get("files", []):
        files.append(
            {
                "id": int(item["id"]),
                "name": str(item["name"]),
                "size": int(item.get("size", 0)),
                "download_url": str(item["download_url"]),
                "md5": str(item.get("computed_md5") or item.get("supplied_md5") or ""),
            }
        )
    return files


def _safe_target(root: Path, name: str) -> Path:
    target = (root / name).resolve()
    root = root.resolve()
    if target != root and root not in target.parents:
        raise ValueError(f"unsafe archive path: {name}")
    return target


def safe_extract(path: str | Path, destination: str | Path) -> list[Path]:
    path = Path(path)
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as archive:
            for member in archive.infolist():
                _safe_target(destination, member.filename)
            archive.extractall(destination)
            extracted = [_safe_target(destination, member.filename) for member in archive.infolist()]
    elif tarfile.is_tarfile(path):
        with tarfile.open(path) as archive:
            for member in archive.getmembers():
                _safe_target(destination, member.name)
                if member.issym() or member.islnk():
                    raise ValueError(f"links are not extracted from QMOF archives: {member.name}")
                if not (member.isfile() or member.isdir()):
                    raise ValueError(f"unsupported special archive member: {member.name}")
            for member in archive.getmembers():
                archive.extract(member, destination)
            extracted = [_safe_target(destination, member.name) for member in archive.getmembers()]
    else:
        raise ValueError(f"unsupported archive format: {path}")
    return extracted


def download_files(
    output_dir: str | Path,
    *,
    names: set[str] | None = None,
    match: str | None = None,
    extract: bool = False,
) -> list[dict]:
    available = list_files()
    regex = re.compile(match) if match else None
    selected = [
        item
        for item in available
        if (names is not None and item["name"] in names)
        or (regex is not None and regex.search(item["name"]))
    ]
    if names is None and regex is None:
        raise ValueError("select files with names or a regular expression")
    if not selected:
        raise ValueError("no QMOF Figshare files matched the selection")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for item in selected:
        target = _safe_target(output_dir, item["name"])
        temporary = target.with_name(target.name + ".part")
        hasher = hashlib.md5()  # nosec B324 - Figshare integrity field, not security use
        request = urllib.request.Request(
            item["download_url"],
            headers={"User-Agent": "qmof-relax-filter/0.1 (+academic software)"},
        )
        with urllib.request.urlopen(request, timeout=120) as response, open(temporary, "wb") as handle:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
                hasher.update(chunk)
        digest = hasher.hexdigest()
        if item["md5"] and digest.lower() != item["md5"].lower():
            temporary.unlink(missing_ok=True)
            raise IOError(f"checksum mismatch for {item['name']}")
        temporary.replace(target)
        extracted_to = None
        if extract and (zipfile.is_zipfile(target) or tarfile.is_tarfile(target)):
            extracted_to = output_dir / target.name.split(".", 1)[0]
            safe_extract(target, extracted_to)
        results.append(
            {
                "name": item["name"],
                "path": str(target),
                "size": target.stat().st_size,
                "md5": digest,
                "extracted_to": str(extracted_to) if extracted_to else None,
            }
        )
    return results
