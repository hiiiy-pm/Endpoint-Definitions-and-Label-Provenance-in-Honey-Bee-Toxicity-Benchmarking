"""Cache a small set of public primary-source provenance records."""
from pathlib import Path
import json
import urllib.request
import urllib.error

OUT=Path(__file__).resolve().parent/'cache'
sources={
 'upstream_ecotox.py':'https://raw.githubusercontent.com/j-adamczyk/ApisTox_dataset/master/dataset_creation/ecotox.py',
 'upstream_processing.py':'https://raw.githubusercontent.com/j-adamczyk/ApisTox_dataset/master/dataset_creation/processing.py',
 'upstream_ppdb_and_bpdb.py':'https://raw.githubusercontent.com/j-adamczyk/ApisTox_dataset/master/dataset_creation/ppdb_and_bpdb.py',
 'upstream_ecotox.csv':'https://zenodo.org/records/13350981/files/ecotox.csv?download=1',
 'upstream_ppdb.csv':'https://zenodo.org/records/13350981/files/ppdb.csv?download=1',
 'oft_zenodo.json':'https://zenodo.org/api/records/19388272',
}
results=[]
for filename,url in sources.items():
    path=OUT/filename
    try:
        if not path.exists():
            with urllib.request.urlopen(url,timeout=30) as r:
                path.write_bytes(r.read())
        results.append({'file':filename,'url':url,'bytes':path.stat().st_size})
    except (urllib.error.URLError,TimeoutError) as e:
        results.append({'file':filename,'url':url,'error':str(e)})
    print(results[-1],flush=True)
(OUT/'provenance_fetch_log.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
