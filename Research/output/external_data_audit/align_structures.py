"""Build auditable structure identities. No source file is modified.

Full standard InChIKey is the primary overlap key. Connectivity and standardized
parent connectivity are conservative sensitivity keys, not identity assertions.
Name/CAS are used only to locate structure records, never as the overlap key.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import time
import urllib.request
import urllib.parse
import urllib.error
from collections import defaultdict

import pandas as pd
from rdkit import Chem, rdBase, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

BASE = (Path(__file__).resolve().parents[2] / 'data/external_raw')
OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
CACHE = OUT/'cache'
NET = CACHE/'pubchem'
NET.mkdir(exist_ok=True)
RDLogger.DisableLog('rdApp.*')

def s(x):
    return '' if pd.isna(x) else str(x).strip()

def norm(x):
    return re.sub(r'[^a-z0-9]', '', s(x).lower())

def safe_key(mol):
    try:
        return Chem.MolToInchiKey(mol) if mol is not None else ''
    except Exception:
        return ''

def structure(smiles='', inchi='', supplied_key=''):
    result = dict(smiles_input=s(smiles), inchi_input=s(inchi), supplied_key=s(supplied_key),
                  canonical_smiles='', inchikey='', connectivity='', parent_connectivity='',
                  fragments=0, structure_status='unresolved', parent_status='unresolved', inchi_smiles_disagree=False)
    mol = Chem.MolFromSmiles(s(smiles)) if s(smiles) else None
    imol = Chem.MolFromInchi(s(inchi)) if s(inchi) else None
    if mol is not None and imol is not None:
        result['inchi_smiles_disagree'] = safe_key(mol) != safe_key(imol)
    if mol is None:
        mol = imol
    if mol is None:
        return result
    key = safe_key(mol)
    if not re.fullmatch('[A-Z]{14}-[A-Z]{10}-[A-Z]', key):
        return result
    parent_key = ''
    try:
        parent = rdMolStandardize.FragmentParent(mol)
        parent = rdMolStandardize.Uncharger().uncharge(parent)
        Chem.RemoveStereochemistry(parent)
        parent_key = safe_key(parent).split('-')[0]
    except Exception:
        pass
    result.update(canonical_smiles=Chem.MolToSmiles(mol), inchikey=key,
                  connectivity=key.split('-')[0],
                  parent_connectivity=parent_key, parent_status='resolved' if parent_key else 'failed',
                  fragments=len(Chem.GetMolFrags(mol)), structure_status='resolved')
    return result

def fetch_pubchem(query, *, offline=False):
    path = NET/(hashlib.sha256(query.encode()).hexdigest()+'.json')
    bundled = CACHE / 'pubchem_records.json'
    if path.exists():
        previous = json.loads(path.read_text(encoding='utf-8'))
    elif bundled.exists():
        previous = json.loads(bundled.read_text(encoding='utf-8')).get(path.name, {})
    else:
        previous = {}
    if previous and '503' not in previous.get('error',''):
        return previous
    if offline:
        raise RuntimeError(f"Offline reconstruction requires a complete, non-retryable bundled PubChem response for {query!r}; no request was sent")
    url = ('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
           + urllib.parse.quote(query, safe='') + '/property/InChIKey,InChI,IsomericSMILES/JSON')
    record = {'query': query, 'url':url, 'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    record['previous_attempts'] = previous.get('previous_attempts',[])
    if previous.get('error'):
        record['previous_attempts'].append(previous['error'])
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=25) as response:
                record['response'] = json.loads(response.read())
            record.pop('error',None)
            break
        except (urllib.error.URLError, TimeoutError) as e:
            record['error'] = str(e)
            if '503' not in str(e):
                break
            if attempt < 2:
                time.sleep(2*(attempt+1))
    path.write_text(json.dumps(record,ensure_ascii=False,indent=2),encoding='utf-8')
    time.sleep(0.25)
    return record

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offline', action='store_true', help='Require existing PubChem responses; never query the network')
    parser.add_argument('--rebuild-structures', action='store_true', help='Recompute OFT reference structures from the extracted source workbook')
    args = parser.parse_args()
    a = pd.read_csv(BASE/'dataset_final.csv')
    a = pd.concat([a,pd.DataFrame([structure(x) for x in a.SMILES])],axis=1)
    a['tier'] = [0 if p==0 else 3 if p==2 else 2 if l==1 else 1 for p,l in zip(a.ppdb_level,a.label)]
    a.to_csv(OUT/'apistox_structures.csv',index=False,encoding='utf-8-sig')
    print('ApisTox',len(a),a.inchikey.nunique(),flush=True)
    path = CACHE/'oft_reference_structures.pkl'
    if path.exists() and not args.rebuild_structures:
        refs = pd.read_pickle(path)
    else:
        refs = pd.read_pickle(CACHE/'REF_SUB.pkl')
        ids = []
        for i, row in refs.iterrows():
            ids.append(structure(row['MolecularStructuralInfo.SmilesNotation'],
                                 row['MolecularStructuralInfo.InChl'],row['MolecularStructuralInfo.InChIKey']))
            if (i+1)%1000 == 0:
                print('OFT structures', i+1,flush=True)
        refs = pd.concat([refs,pd.DataFrame(ids)],axis=1)
        refs.to_pickle(path)
    print('OFT', len(refs),refs.structure_status.value_counts().to_dict(),flush=True)
    # Only unambiguous names/CAS link to a reference structure. Retain ambiguity.
    byname, bycas = defaultdict(list), defaultdict(list)
    namecols = ['ReferenceSubstanceName','IupacName','CAS name','COM NAME [EFSA OFT2.0]',
                'Name','PARAM NAME','SUB NAME [EFSA OFT2.0]']
    for i,row in refs.iterrows():
        if not row['inchikey']:
            continue
        for col in namecols:
            for name in s(row[col]).split('|'):
                if name:
                    byname[norm(name)].append(i)
        for col in ['Inventory.CASNumber','CAS number']:
            for cas in s(row[col]).split('|'):
                if cas:
                    bycas[norm(cas)].append(i)
    key_fields = list(structure().keys())
    api_keysets = {k:set(a[k])-{''} for k in ['inchikey','connectivity','parent_connectivity']}

    def resolve(query, mode='name'):
        indices = (byname if mode=='name' else bycas).get(norm(query),[])
        hit = refs.loc[sorted(set(indices))]
        if len(hit) and hit.inchikey.nunique() == 1:
            chosen = hit.iloc[0]
            return {**{k:chosen[k] for k in key_fields},'mapping_method':'OFT exact '+mode,
                    'mapping_reference':'|'.join(hit['Document UUID']), 'mapping_query':query,
                    'mapping_CAS':s(chosen['Inventory.CASNumber'])}
        fetched = fetch_pubchem(query, offline=args.offline)
        props = fetched.get('response',{}).get('PropertyTable',{}).get('Properties',[])
        if props and len({x.get('InChIKey') for x in props})==1:
            p = props[0]
            z = structure(p.get('SMILES',p.get('IsomericSMILES','')),p.get('InChI',''),p.get('InChIKey',''))
            return {**z,'mapping_method':'PubChem '+mode, 'mapping_reference':fetched['url'],
                    'mapping_query':query, 'pubchem_cid':p.get('CID'),
                    'mapping_CAS':query if mode=='CAS' else ''}
        return {**structure(),'mapping_method':'unresolved/ambiguous','mapping_reference':fetched['url'],
                'mapping_query':query,'candidate_records':len(props),'local_candidate_keys':hit.inchikey.nunique(),
                'candidate_connectivities':'|'.join(sorted({x.get('InChIKey','').split('-')[0] for x in props})),
                'mapping_CAS':''}

    p = pd.read_pickle(CACHE/'plos_s1.pkl')
    p = p[p.Pesticide.notna()].copy()
    mapped=[]
    for i,row in p.iterrows():
        mapped.append(resolve(row.Pesticide))
        if (i+1)%25==0:
            print('PLOS identities',i+1,flush=True)
    p=pd.concat([p.reset_index(drop=True),pd.DataFrame(mapped)],axis=1)
    for key, values in api_keysets.items():
        p['overlap_'+key]=p[key].isin(values) & p[key].ne('')
    p.to_csv(OUT/'plos_identity_alignment.csv',index=False,encoding='utf-8-sig')
    print('PLOS',p.structure_status.value_counts().to_dict(),p.mapping_method.value_counts().to_dict(),flush=True)
    u = pd.read_csv(OUT/'excluded_all_parsed.csv')
    u=u[u['Exclusion reason'].eq('Unspecifed toxicity level')].drop_duplicates('CAS')[['CAS','name']]
    mapped=[]
    for _,row in u.iterrows():
        mapped.append(resolve(row.CAS, 'CAS'))
    u=pd.concat([u.reset_index(drop=True),pd.DataFrame(mapped)],axis=1)
    for key,values in api_keysets.items():
        u['overlap_'+key]=u[key].isin(values) & u[key].ne('')
    u.to_csv(OUT/'excluded_uncertain_identity_alignment.csv',index=False,encoding='utf-8-sig')
    print('UNCERTAIN',u.structure_status.value_counts().to_dict(),flush=True)
    meta={'rdkit':rdBase.rdkitVersion,'pandas':pd.__version__,
          'download_matches_frozen_table':pd.read_csv(BASE/'dataset_final.csv').equals(pd.read_csv(ROOT/'data/raw/apistox.csv')),
          'primary_key':'standard full InChIKey computed from supplied SMILES, falling back to supplied InChI',
          'parent_sensitivity':'RDKit FragmentParent, Uncharger, remove stereochemistry; first InChIKey block',
          'unresolved_policy':'Excluded from identity-overlap denominators; never called novel'}
    assert meta['download_matches_frozen_table'], 'External benchmark copy differs from the frozen raw benchmark'
    (OUT/'structure_methods.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')

if __name__=='__main__':
    main()
