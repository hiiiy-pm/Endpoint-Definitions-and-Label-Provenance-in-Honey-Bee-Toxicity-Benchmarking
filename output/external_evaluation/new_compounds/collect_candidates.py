"""Read-only source collection and structural exclusion against frozen ApisTox.

No model fitting or prediction. Candidate list fixed before viewing predictions.
"""
from pathlib import Path
import sys, json, urllib.request, urllib.parse, hashlib, concurrent.futures
import pandas as pd

OUT=Path(__file__).resolve().parent
AUDIT=OUT.parents[1]/'external_data_audit'
OLD=OUT.parents[1]/'additional_external_sources'
sys.path.insert(0,str(AUDIT))
from align_structures import structure
A=pd.read_csv(AUDIT/'apistox_structures.csv').fillna('')
CANDIDATES=['cyclobutrifluram','isocycloseram','dimpropyridaz','acynonapyr','oxazosulfyl','fluhexafon','inpyrfluxam','spiropidion']
SOURCES={
 'acynonapyr_pmra_2026.pdf':'https://publications.gc.ca/collections/collection_2026/sc-hc/h113-9/H113-9-2026-3-eng.pdf',
 'spiropidion_pmra_2026.pdf':'https://publications.gc.ca/collections/collection_2026/sc-hc/h113-9/H113-9-2026-9-eng.pdf',
 'cyclobutrifluram_epa_2025.pdf':'https://downloads.regulations.gov/EPA-HQ-OPP-2022-0003-0026/content.pdf',
 'inpyrfluxam_epa.pdf':'https://downloads.regulations.gov/EPA-HQ-OPP-2018-0038-0025/content.pdf',
 'inpyrfluxam_developer_2020.pdf':'https://www.sumitomo-chem.co.jp/rd/report/files/docs/2020J_1.pdf',
 'oxazosulfyl_developer_2021.pdf':'https://www.sumitomo-chem.co.jp/english/rd/report/2021E_2.pdf',
 'isocycloseram_who_2025.pdf':'https://extranet.who.int/prequal/sites/default/files/doc_parts/WHOVC-SP_Isocycloseram_2025.1.pdf',
}

def identity(query):
    url='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.parse.quote(query)+'/property/InChIKey,InChI,IsomericSMILES/JSON'
    cache=OUT/(query+'_pubchem.json')
    old=OLD/(query+'_pubchem.json')
    data={}
    try:
        if cache.exists(): data=json.loads(cache.read_text(encoding='utf-8'))
        elif old.exists(): data=json.loads(old.read_text(encoding='utf-8'))
        else:
            with urllib.request.urlopen(url,timeout=25) as resp:data={'query_url':url,'response':json.loads(resp.read())}
    except Exception as e:data={'query_url':url,'error':str(e)}
    cache.write_text(json.dumps(data,indent=2),encoding='utf-8')
    props=data.get('response',data).get('PropertyTable',{}).get('Properties',[])
    rows=[]
    for p in props:
        z=structure(p.get('SMILES',p.get('IsomericSMILES','')),p.get('InChI',''),p.get('InChIKey',''))
        r={'candidate':query,'PubChem_CID':p.get('CID'),'structure_source':url,**z}
        for k in ['inchikey','connectivity','parent_connectivity']:
            m=A[A[k].eq(z[k])] if z[k] else A.iloc[:0]
            r['overlap_'+k]=bool(len(m));r['matched_names_'+k]='|'.join(m.name)
        rows.append(r)
    return rows or [{'candidate':query,'structure_status':'unresolved','error':data.get('error','No properties')}]

def source(item):
    name,url=item;p=OUT/'sources'/name;p.parent.mkdir(exist_ok=True)
    row={'file':str(p),'url':url}
    try:
        if not p.exists():
            with urllib.request.urlopen(url,timeout=30) as r:b=r.read()
            if not b.startswith(b'%PDF'):raise ValueError('response is not PDF')
            p.write_bytes(b)
        row.update(status='downloaded',sha256=hashlib.sha256(p.read_bytes()).hexdigest())
        try:
            from pypdf import PdfReader
            pdf=PdfReader(str(p));t='\n'.join(f'\n=== PDF PAGE {i+1} ===\n'+(x.extract_text() or '') for i,x in enumerate(pdf.pages))
            p.with_suffix('.txt').write_text(t,encoding='utf-8');row['pages']=len(pdf.pages)
        except ImportError:row['text_extraction']='pypdf unavailable'
    except Exception as e:row.update(status='unavailable',error=str(e))
    return row

if __name__=='__main__':
    rows=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for r in pool.map(identity,CANDIDATES):rows.extend(r)
    df=pd.DataFrame(rows);df.to_csv(OUT/'candidate_structure_exclusion.csv',index=False,encoding='utf-8-sig')
    print(df[['candidate','PubChem_CID','inchikey','overlap_parent_connectivity']].to_string(index=False),flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:sr=list(pool.map(source,SOURCES.items()))
    pd.DataFrame(sr).to_csv(OUT/'source_download_manifest.csv',index=False,encoding='utf-8-sig')
    print(pd.DataFrame(sr)[['file','status']].to_string(index=False))
