"""Current PubChem identity checks; never changes the frozen benchmark.

Raw responses are timestamped current evidence, not reconstructed 2024 snapshots.
Run with an environment containing requests and, for compare, RDKit.
"""
from pathlib import Path
import csv, json, hashlib, datetime, time, argparse, collections
from urllib.parse import quote

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
RAW = RESEARCH / 'data/raw/apistox.csv'
COMMIT = '57c2c05de293a0f54c119744eb3ca9cb03bc030d'
PROPERTIES = 'InChIKey,IsomericSMILES,MolecularFormula,IUPACName'
PUG = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'
ROWS = list(csv.DictReader(RAW.open(encoding='utf-8-sig')))
for i, row in enumerate(ROWS): row['Index'] = i

def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding='utf-8')

def fetch(url, stem, params=None, force=False):
    import requests
    out = HERE / 'responses' / (stem + '.body')
    meta = HERE / 'responses' / (stem + '.json')
    if out.exists() and meta.exists() and not force:
        m = json.loads(meta.read_text(encoding='utf-8'))
        if m.get('status') == 200:
            try: return json.loads(out.read_text(encoding='utf-8'))
            except json.JSONDecodeError: return None
    for attempt in range(3):
        started = stamp()
        try:
            response = requests.get(url, params=params, timeout=30)
            content = response.content
            out.write_bytes(content)
            write_json(meta, {'request_started_utc': started, 'response_completed_utc': stamp(),
                'url': response.url, 'status': response.status_code,
                'response_headers': dict(response.headers), 'sha256': hashlib.sha256(content).hexdigest(),
                'bytes': len(content), 'attempt': attempt + 1,
                'scope': 'current_public_response_not_a_2024_snapshot'})
            time.sleep(.25)
            if response.status_code == 200:
                try: return response.json()
                except ValueError: return None
            if response.status_code not in (429, 500, 502, 503, 504): return None
        except requests.RequestException as error:
            write_json(meta, {'request_started_utc': started, 'response_completed_utc': stamp(),
                'url': url, 'status': 'request_failed', 'error': str(error), 'attempt': attempt+1,
                'scope': 'current_public_response_not_a_2024_snapshot'})
        time.sleep(attempt + 1)
    return None

def fetch_all():
    cids = sorted({int(row['CID']) for row in ROWS})
    for i in range(0, len(cids), 40):
        batch = cids[i:i+40]
        fetch(PUG + 'cid/' + ','.join(map(str,batch)) + '/property/' + PROPERTIES + '/JSON',
              'cid_batch_%04d' % i)
        print('CID batch', i, 'of', len(cids), flush=True)
    for filename in ['dataset_creation/pubchem.py', 'dataset_creation/ecotox.py',
                     'dataset_creation/ppdb_and_bpdb.py', 'create_dataset.py', 'create_dataset_splits.py']:
        url = 'https://raw.githubusercontent.com/j-adamczyk/ApisTox_dataset/' + COMMIT + '/' + filename
        stem = 'upstream_' + filename.replace('/', '_')
        fetch(url, stem)
        src = HERE / 'responses' / (stem + '.body')
        if src.exists(): (HERE / 'upstream' / filename.replace('/', '_')).write_bytes(src.read_bytes())

def fetch_known():
    grouped = collections.defaultdict(list)
    for row in ROWS: grouped[row['CID']].append(row)
    known = [row for group in grouped.values() if len(group)>1 for row in group]
    for row in known:
        fetch(PUG + 'name/' + quote(row['CAS'], safe='') + '/property/' + PROPERTIES + '/JSON',
              'cas_' + row['CAS'])
        print('CAS candidate', row['Index'], row['CAS'], flush=True)
    candidate_cids = {int(row['CID']) for row in known}
    for row in known:
        p=HERE/'responses'/('cas_'+row['CAS']+'.body')
        if p.exists():
            try: candidate_cids.update(x['CID'] for x in json.loads(p.read_text(encoding='utf-8'))['PropertyTable']['Properties'])
            except (KeyError, ValueError): pass
    for cid in sorted(candidate_cids):
        query={'download':'*','collection':'literature','where':{'ands':[{'cid':str(cid)}]},
               'order':['articlepubdate,asc'],'start':0,'limit':5,
               'downloadfilename':f'pubchem_cid_{cid}_literature','nullatbottom':1}
        fetch('https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi', 'literature_'+str(cid),
              {'infmt':'json','outfmt':'json','query':json.dumps(query)})
        fetch('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/'+str(cid)+'/JSON/',
              'pugview_'+str(cid))
        print('CID provenance', cid, flush=True)

def fetch_candidates():
    flagged=list(csv.DictReader((HERE/'current_identifier_flagged.csv').open(encoding='utf-8-sig')))
    for row in flagged:
        cas=row['CAS'].strip()
        fetch(PUG+'name/'+quote(cas,safe='')+'/property/'+PROPERTIES+'/JSON','cas_'+cas)
        print('Flagged CAS candidate',row['Index'],cas,flush=True)
    fetch('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/substance/481109147/JSON/',
          'metiram_reference_substance_481109147')
    fetch('https://www.epa.gov/sites/default/files/documents/2001ebdcnew.pdf','epa_metiram_zineb_2001')

def write_csv(name, rows):
    if not rows: return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (HERE/name).open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)

def compare():
    from rdkit import Chem, rdBase
    from rdkit.Chem.MolStandardize import rdMolStandardize
    props={}
    for p in sorted((HERE/'responses').glob('cid_batch_*.body')):
        try:
            for r in json.loads(p.read_text(encoding='utf-8'))['PropertyTable']['Properties']: props[r['CID']]=r
        except (ValueError,KeyError): pass
    def keys(smiles):
        mol=Chem.MolFromSmiles(smiles)
        if mol is None:return ('','','')
        full=Chem.MolToInchiKey(mol)
        try: parent=Chem.MolToInchiKey(rdMolStandardize.FragmentParent(mol))
        except Exception: parent=''
        return full,full.split('-')[0],parent
    ledger=[]
    for row in ROWS:
        key,conn,parent=keys(row['SMILES']); current=props.get(int(row['CID']),{})
        ckey=current.get('InChIKey',''); csmiles=current.get('SMILES',current.get('IsomericSMILES',''))
        computed,cconn,cparent=keys(csmiles) if csmiles else ('','','')
        if not current: status='current_CID_response_unavailable'
        elif key==ckey: status='full_InChIKey_match'
        elif conn==ckey.split('-')[0]: status='connectivity_match_full_key_difference'
        elif parent and parent==cparent: status='fragment_parent_match_full_structure_difference'
        elif parent and cparent and parent.split('-')[0]==cparent.split('-')[0]: status='fragment_parent_connectivity_match'
        else: status='structure_mismatch_requires_adjudication'
        ledger.append(dict(row,benchmark_inchikey=key,benchmark_parent_inchikey=parent,
            current_cid_inchikey=ckey,current_cid_smiles=csmiles,current_cid_iupac=current.get('IUPACName',''),
            current_parent_inchikey=cparent,comparison=status,
            current_response_internal_key_matches=str(computed==ckey),
            evidence_epoch='current_2026_response;not_historical_snapshot',
            frozen_fields_changed=False))
    write_csv('cid_structure_ledger_1035.csv',ledger)
    write_csv('cid_structure_review_candidates.csv',[r for r in ledger if r['comparison']!='full_InChIKey_match'])
    summary={'created_utc':stamp(),'rdkit_version':rdBase.rdkitVersion,'benchmark_rows':len(ROWS),
        'unique_CIDs':len({r['CID'] for r in ROWS}),'current_CID_property_records':len(props),
        'status_counts':dict(collections.Counter(r['comparison'] for r in ledger)),
        'source_sha256':hashlib.sha256(RAW.read_bytes()).hexdigest(),
        'claim_limit':'Current identifier consistency check; differences are candidates for review, not automatically corrected benchmark structures, labels, dates or use flags.'}
    write_json(HERE/'identity_summary.json',summary)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['fetch_all','fetch_known','fetch_candidates','compare'])
    globals()[parser.parse_args().mode]()
