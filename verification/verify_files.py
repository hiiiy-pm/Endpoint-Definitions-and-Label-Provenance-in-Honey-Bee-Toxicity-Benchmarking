"""Verify the files shipped in this repository without fitting models."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
manifest = json.loads((ROOT / "FILE_MANIFEST.json").read_text(encoding="utf-8"))
missing, changed = [], []
for entry in manifest["files"]:
    path = ROOT / entry["path"]
    if not path.is_file():
        missing.append(entry["path"])
    elif hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
        changed.append(entry["path"])
report = {
    "files_checked": len(manifest["files"]),
    "missing": missing,
    "changed": changed,
    "status": "PASS" if not missing and not changed else "FAIL",
}
print(json.dumps(report, ensure_ascii=False, indent=2))
if missing or changed:
    raise SystemExit(1)
