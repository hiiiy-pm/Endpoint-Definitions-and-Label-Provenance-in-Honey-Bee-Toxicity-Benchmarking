#!/usr/bin/env python3
"""Build and verify the public research archive deterministically."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import sys
import zipfile

TAG = "research-reproducibility-2026-10-24"
SOURCE_DIRECTORIES = ("Research",)
ASSET_NAMES = ("Research.zip",)

ZIP_TIMESTAMP = (2026, 10, 24, 0, 0, 0)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path, name: str, key: str) -> dict:
    return {key: name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}


def source_files(repository: Path) -> dict[str, Path]:
    result = {}
    for directory in SOURCE_DIRECTORIES:
        root = repository / directory
        if not root.is_dir() or root.is_symlink():
            raise ValueError(f"Required source directory missing or linked: {directory}")
        for path in root.rglob("*"):
            if path.is_symlink():
                raise ValueError(f"Symbolic links are not permitted: {path.relative_to(repository)}")
            if path.is_file():
                result[path.relative_to(repository).as_posix()] = path
    return dict(sorted(result.items()))


def check_record(record: dict, path_key: str) -> None:
    if not isinstance(record, dict) or set(record) != {path_key, "bytes", "sha256"}:
        raise ValueError("Invalid manifest file record")
    name = record[path_key]
    if not isinstance(name, str) or "\\" in name or "\x00" in name:
        raise ValueError("Invalid manifest path")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != name:
        raise ValueError(f"Unsafe or noncanonical manifest path: {name}")
    if not isinstance(record["bytes"], int) or record["bytes"] < 0:
        raise ValueError(f"Invalid byte count: {name}")
    digest = record["sha256"]
    if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
        raise ValueError(f"Invalid SHA-256 digest: {name}")


def read_manifest(path: Path) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or manifest.get("tag") != TAG:
        raise ValueError("Unexpected manifest schema or release tag")
    if not isinstance(manifest.get("files"), list) or not isinstance(manifest.get("assets"), list):
        raise ValueError("Manifest must contain source files and asset records")
    for record in manifest["files"]:
        check_record(record, "path")
        if PurePosixPath(record["path"]).parts[0] not in SOURCE_DIRECTORIES:
            raise ValueError(f"Unexpected source directory: {record['path']}")
    paths = [record["path"] for record in manifest["files"]]
    if len(set(paths)) != len(paths) or manifest.get("source_member_count") != len(paths):
        raise ValueError("Manifest source count or unique paths do not match")
    for record in manifest["assets"]:
        check_record(record, "name")
    names = [record["name"] for record in manifest["assets"]]
    if len(names) != len(ASSET_NAMES) or set(names) != set(ASSET_NAMES):
        raise ValueError("Manifest must contain exactly the expected research archive")
    return manifest


def verify_sources(repository: Path, manifest: dict) -> dict[str, Path]:
    actual = source_files(repository)
    expected = {record["path"]: record for record in manifest["files"]}
    if set(actual) != set(expected):
        missing = sorted(set(expected) - set(actual))
        extra = sorted(set(actual) - set(expected))
        raise ValueError(f"Source file set differs from manifest; missing={missing}, extra={extra}")
    for name, path in actual.items():
        if file_record(path, name, "path") != expected[name]:
            raise ValueError(f"Source bytes or SHA-256 differ from manifest: {name}")
    return actual


def write_zip(path: Path, members: list[tuple[str, Path]]) -> None:
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6, allowZip64=True) as archive:
        for name, source in sorted(members, key=lambda item: item[0]):
            info = zipfile.ZipInfo(name, date_time=ZIP_TIMESTAMP)
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, source.read_bytes(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=6)
    with zipfile.ZipFile(path) as archive:
        bad_member = archive.testzip()
        if bad_member:
            raise ValueError(f"ZIP integrity check failed: {path.name}: {bad_member}")


def build_assets(repository: Path, output: Path, files: dict[str, Path]) -> None:
    output = output.resolve()
    for directory in SOURCE_DIRECTORIES:
        if output.is_relative_to((repository / directory).resolve()):
            raise ValueError("Output directory cannot be within a research source directory")
    output.mkdir(parents=True, exist_ok=True)
    write_zip(output / "Research.zip", [(name, path) for name, path in files.items()])


def verify_assets(output: Path, manifest: dict) -> None:
    actual_names = {path.name for path in output.iterdir() if path.is_file()}
    if actual_names != set(ASSET_NAMES):
        raise ValueError("Output directory does not contain exactly the research archive")
    for record in manifest["assets"]:
        path = output / record["name"]
        if path.is_symlink() or not path.is_file() or file_record(path, record["name"], "name") != record:
            raise ValueError(f"Release asset differs from expected bytes or SHA-256: {record['name']}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prepare-metadata", action="store_true", help="Generate the source and asset manifest locally")
    args = parser.parse_args()
    repository = args.repository.resolve()
    manifest_path = args.manifest or repository / ".github/release/research_metadata.json"
    output = args.output_dir.resolve()
    if args.prepare_metadata:
        files = source_files(repository)
        manifest = {
            "schema_version": 1,
            "tag": TAG,
            "source_member_count": len(files),
            "files": [file_record(path, name, "path") for name, path in files.items()],
        }
        build_assets(repository, output, files)
        manifest["assets"] = [file_record(output / name, name, "name") for name in ASSET_NAMES]
        verify_assets(output, manifest)
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    else:
        manifest = read_manifest(manifest_path)
        files = verify_sources(repository, manifest)
        build_assets(repository, output, files)
        verify_assets(output, manifest)
    print(f"Verified {manifest['source_member_count']} source files and {len(ASSET_NAMES)} release assets.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, zipfile.BadZipFile) as exc:
        print(f"Release asset validation failed: {exc}", file=sys.stderr)
        raise SystemExit(1)
