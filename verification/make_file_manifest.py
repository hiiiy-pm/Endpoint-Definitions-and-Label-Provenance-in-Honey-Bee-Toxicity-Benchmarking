"""Write FILE_MANIFEST.json and SHA256SUMS.txt for the files of this repository.

The file list is taken from git (tracked files plus new files that are not
ignored), so build products and local template assets are excluded. The
checksum file uses LF line endings so that ``sha256sum -c SHA256SUMS.txt`` works
on Linux and macOS as well as with verification/verify_files.py.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
EXCLUDE = {'FILE_MANIFEST.json', 'SHA256SUMS.txt'}

listed = subprocess.run(['git', 'ls-files', '--cached', '--others', '--exclude-standard', '-z'],
                        cwd=ROOT, check=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
paths = sorted({p for p in listed.split('\0') if p and p not in EXCLUDE and (ROOT / p).is_file()})
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
