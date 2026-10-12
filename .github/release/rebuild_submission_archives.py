"""Rebuild submitted ZIPs only when all source and final archive hashes match.

Standard library only. No downloads, model fitting, or GitHub writes.
Publication callers MUST require this script's successful exit before uploads.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import os
import platform
import re
import sys
import zipfile
import zlib


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def source_path(root, name):
    archive_path = PurePosixPath(name)
    require(len(archive_path.parts) > 1 and archive_path.parts[0] == "Research" and not archive_path.is_absolute(), "Unexpected ZIP member prefix")
    require(".." not in archive_path.parts and "\\" not in name, "Unsafe ZIP member path")
    path = root.joinpath(*archive_path.parts[1:]).resolve()
    require(path.is_relative_to(root) and path.is_file(), f"Source file missing or outside root: {name}")
    return path


def verify_research_manifest(root, metadata):
    """Every archived research file must be covered, except the two self-indexes.

    Repository-only .github/release helpers and workflow files are not research
    members and therefore do not enter either the manifest or the four ZIPs.
    """
    entries = [entry for archive in metadata["archives"] for entry in archive["entries"]]
    archive_files = {entry["name"].removeprefix("Research/"): entry for entry in entries}
    require(len(archive_files) == metadata["source_member_count"], "Metadata member count differs")
    manifest = json.loads(source_path(root, "Research/FILE_MANIFEST.json").read_text(encoding="utf-8-sig"))
    rows = manifest["files"]
    by_path = {row["path"]: row for row in rows}
    require(len(by_path) == len(rows) == manifest["file_count"] == metadata["research_manifest_entry_count"], "Manifest count or duplicate path mismatch")
    exclusions = {"FILE_MANIFEST.json", "SHA256SUMS.txt"}
    require(set(manifest["manifest_excludes"]) == exclusions, "Unexpected research manifest exclusions")
    require(set(by_path) == set(archive_files) - exclusions, "Research manifest and archive membership differ")
    require(exclusions.issubset(archive_files), "Required research index files are absent")
    for name, row in by_path.items():
        archived = archive_files[name]
        require(row["bytes"] == archived["file_size"] and row["sha256"] == archived["sha256"], f"Manifest and archive metadata disagree: {name}")
        path = source_path(root, "Research/" + name)
        require(path.stat().st_size == row["bytes"] and sha256(path.read_bytes()) == row["sha256"], f"Research manifest verification failed: {name}")
    checksums = {}
    for line in source_path(root, "Research/SHA256SUMS.txt").read_text(encoding="utf-8-sig").splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        require(match is not None, "Malformed research SHA256SUMS line")
        digest, name = match.groups()
        require(name not in checksums, "Duplicate research SHA256SUMS path")
        checksums[name] = digest
    require(checksums == {name: row["sha256"] for name, row in by_path.items()}, "Research SHA256SUMS and manifest differ")
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--research-root", type=Path, required=True, help="Repository root containing README.md, code/, data/, etc.")
    parser.add_argument("--metadata", type=Path, default=Path(__file__).with_name("submission_zip_metadata.json"))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.research_root.resolve(strict=True)
    out = args.output_dir.resolve()
    require(out != root and not root.is_relative_to(out), "Output directory must not contain or equal the source root")
    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    require(metadata["schema_version"] == 2 and len(metadata["archives"]) == 4, "Unrecognized archive metadata")
    require(metadata["compression_level"] == 9, "Unexpected compression level")
    require(isinstance(metadata["source_member_count"], int) and metadata["source_member_count"] > 2, "Invalid explicit source member count")
    names = [entry["name"] for archive in metadata["archives"] for entry in archive["entries"]]
    require(len(names) == len(set(names)) == metadata["source_member_count"], "Archive partition is not complete and disjoint")
    require(len({a["name"] for a in metadata["archives"]}) == 4, "Duplicate archive names")
    manifest_count = verify_research_manifest(root, metadata)
    # Check every source before creating output. Ignore repository-only release helpers.
    for archive in metadata["archives"]:
        for entry in archive["entries"]:
            path = source_path(root, entry["name"])
            require(path.stat().st_size == entry["file_size"], f"Source size mismatch: {entry['name']}")
            require(sha256(path.read_bytes()) == entry["sha256"], f"Source SHA256 mismatch: {entry['name']}")
    out.mkdir(parents=True, exist_ok=True)
    runtime = {"python": platform.python_version(), "zlib_compiled": zlib.ZLIB_VERSION, "zlib_runtime": zlib.ZLIB_RUNTIME_VERSION, "platform": sys.platform}
    print(json.dumps({"runtime": runtime, "recorded_runtime": metadata["captured_runtime"], "source_members_verified": len(names), "research_manifest_entries_verified": manifest_count}), flush=True)
    temporary = []
    verified = []
    try:
        for archive in metadata["archives"]:
            name = archive["name"]
            require(Path(name).name == name and name.startswith("Research_") and name.endswith(".zip"), "Unsafe archive output name")
            temp = out / (name + ".pending")
            require(not temp.exists(), f"Refusing to replace pending output: {temp}")
            temporary.append(temp)
            with zipfile.ZipFile(temp, "w", zipfile.ZIP_DEFLATED, compresslevel=9, allowZip64=True) as result:
                result.comment = bytes.fromhex(archive["archive_comment_hex"])
                for entry in archive["entries"]:
                    info = zipfile.ZipInfo(entry["name"], tuple(entry["date_time"]))
                    for field in ("create_system", "create_version", "extract_version", "reserved", "flag_bits", "volume", "internal_attr", "external_attr", "compress_type"):
                        setattr(info, field, entry[field])
                    info.comment = bytes.fromhex(entry["comment_hex"])
                    info.extra = bytes.fromhex(entry["extra_hex"])
                    info.file_size = entry["file_size"]
                    result.writestr(info, source_path(root, entry["name"]).read_bytes(), compress_type=entry["compress_type"], compresslevel=9)
            actual = sha256(temp.read_bytes())
            require(temp.stat().st_size == archive["bytes"] and actual == archive["sha256"], f"ARCHIVE MISMATCH; DO NOT PUBLISH: {name}; expected={archive['sha256']}; actual={actual}; runtime={runtime}")
            with zipfile.ZipFile(temp) as check:
                require(check.testzip() is None, f"ZIP CRC mismatch: {name}")
                require(check.namelist() == [e["name"] for e in archive["entries"]], f"ZIP ordering mismatch: {name}")
            verified.append({"file": name, "bytes": temp.stat().st_size, "sha256": actual, "members": len(archive["entries"]), "byte_identical_to_submission": True})
            print(json.dumps(verified[-1]), flush=True)
        # No .zip publication candidates exist until every archive has passed.
        for temp, archive in zip(temporary, metadata["archives"], strict=True):
            os.replace(temp, out / archive["name"])
        report = {"status": "PASS", "runtime": runtime, "source_members_verified": len(names), "research_manifest_entries_verified": manifest_count, "metadata_sha256": sha256(args.metadata.read_bytes()), "archives": verified, "experiments_run": False, "github_operations": False}
        (out / "archive_rebuild_verification.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        (out / "RESEARCH_ARCHIVE_SHA256SUMS.txt").write_text("".join(f"{r['sha256']}  {r['file']}\n" for r in verified), encoding="utf-8", newline="\n")
        print("PASS: all four archives are byte-identical to the submitted files.", flush=True)
    finally:
        for temp in temporary:
            if temp.exists():
                temp.unlink()


if __name__ == "__main__":
    main()
