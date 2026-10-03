"""Rebuild the delivered-file manifest after verification and artifact generation.

Use git's ignore rules for generated scratch material; no commit is created.
Run only in an authored checkout, not merely to hide unexpected file changes.
"""
from pathlib import Path
import hashlib
import json
import subprocess

ROOT=Path(__file__).resolve().parents[1]
raw=subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=ROOT)
paths=sorted(set(raw.decode('utf-8').split('\0'))-{'','FILE_MANIFEST.json','SHA256SUMS.txt'})
rows=[]
for rel in paths:
    p=ROOT/rel
    if not p.is_file():
        raise FileNotFoundError(rel)
    b=p.read_bytes()
    rows.append({'path':rel,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
(ROOT/'FILE_MANIFEST.json').write_text(json.dumps({'files':rows},indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
(ROOT/'SHA256SUMS.txt').write_text(''.join(f"{x['sha256']}  {x['path']}\n" for x in rows),encoding='utf-8')
print('Manifest entries:',len(rows))
