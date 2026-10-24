"""Project a frozen curator table to permitted identifiers and predict once.

This runner records hashes/timestamps before inference, never reads toxicity
labels, and prints only completion counts, not prediction outcomes.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import sys
import pandas as pd

OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(OUT/'model_restore'))
from predict_external import predict_structure_table

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def utc():
    return datetime.now(timezone.utc).isoformat()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',required=True)
    ap.add_argument('--source-manifest',required=True)
    ap.add_argument('--batch',required=True)
    args=ap.parse_args()
    if not args.batch.replace('_','').isalnum():
        raise ValueError('Batch must be a simple alphanumeric identifier')
    target=OUT/'predictions'/args.batch
    target.mkdir(parents=True,exist_ok=True)
    log_path=target/'prediction_run_manifest.json'
    if log_path.exists():
        old=json.loads(log_path.read_text(encoding='utf-8'))
        if old.get('completed'):
            raise RuntimeError('Completed prediction batch is immutable; verify or use a new version')
    source=Path(args.input).resolve()
    source_manifest=Path(args.source_manifest).resolve()
    meta=json.loads(source_manifest.read_text(encoding='utf-8'))
    expected=meta.get('sha256')
    if not expected:
        raise ValueError('Curator manifest must explicitly contain sha256 for input structure table')
    if sha(source)!=expected:
        raise ValueError('Frozen structure table checksum mismatch')
    table=pd.read_csv(source,keep_default_na=False)
    forbidden=[c for c in table if any(x in c.lower() for x in ['ld50','toxicity','endpoint','label_le','y_1','y_11','y_100','mortality'])]
    if forbidden:
        raise ValueError(f'Frozen structure table contains endpoint fields: {forbidden}')
    allowed=[c for c in ['compound_id','SMILES','name','inchikey','CAS','source'] if c in table]
    if not {'compound_id','SMILES'}.issubset(allowed):
        raise ValueError('Missing required structure fields')
    projection=table[allowed]
    inference_input=target/'structure_only_input.csv'
    projection.to_csv(inference_input,index=False,encoding='utf-8-sig')
    model_manifest=OUT/'model_restore/restoration_manifest.json'
    record={'started_utc':utc(),'completed':False,'source_path':str(source),
        'source_sha256':sha(source),'source_manifest_sha256':sha(source_manifest),
        'source_last_modified_utc':datetime.fromtimestamp(source.stat().st_mtime,timezone.utc).isoformat(),
        'structure_only_input_sha256':sha(inference_input),
        'restoration_manifest_sha256':sha(model_manifest),'runner_sha256':sha(__file__),
        'n_structures':len(projection),'fit_calls_in_runner':0,'external_labels_read':False,
        'require_new_structures':True,'columns_passed_to_inference':allowed}
    log_path.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    predictions=predict_structure_table(projection,model_dir=OUT/'model_restore',require_new_structures=True)
    path=target/'new_compound_predictions.csv'
    predictions.to_csv(path,index=False,encoding='utf-8-sig')
    record.update(completed=True,completed_utc=utc(),n_prediction_rows=len(predictions),
                  prediction_sha256=sha(path),source_unchanged_after_inference=sha(source)==expected)
    assert record['source_unchanged_after_inference']
    log_path.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    print(f'Completed frozen batch {args.batch}: {len(projection)} structures, {len(predictions)} predictions. No fitting or labels in inference.')

if __name__=='__main__':
    main()
