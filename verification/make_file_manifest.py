"""Refresh the authored archive manifest without requiring a Git checkout.

Run after intentional edits and verification, not to conceal unexpected changes.
Historical run manifests are evidence of those runs and are not rewritten here.
"""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
EXCLUDE = {'FILE_MANIFEST.json', 'SHA256SUMS.txt'}

SKIP_DIRS = {'.git', '.venv', 'venv', '__pycache__', '.pytest_cache'}
paths = sorted(p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
               if p.is_file() and p.relative_to(ROOT).as_posix() not in EXCLUDE
               and not any(part in SKIP_DIRS for part in p.relative_to(ROOT).parts)
               and p.suffix not in {'.pyc', '.pyo'})
files = []
for rel in paths:
    data = (ROOT / rel).read_bytes()
    files.append({'path': rel, 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()})
manifest = {'files': files, 'file_count': len(files),
            'manifest_excludes': sorted(EXCLUDE)}
(ROOT / 'FILE_MANIFEST.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + '\n',
                                           encoding='utf-8', newline='\n')
with open(ROOT / 'SHA256SUMS.txt', 'w', encoding='utf-8', newline='\n') as handle:
    for entry in files:
        handle.write(f"{entry['sha256']}  {entry['path']}\n")
print(f'{len(files)} files listed')
