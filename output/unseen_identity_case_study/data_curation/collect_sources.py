"""Bounded candidate collection. No training or prediction files are accessed."""
from pathlib import Path
import sys,json,hashlib,urllib.request,urllib.parse,concurrent.futures
import pandas as pd
OUT=Path(__file__).resolve().parent
OUTPUT=OUT.parents[1]
OLD=OUTPUT/'external_evaluation/new_compounds'
AUDIT=OUTPUT/'external_data_audit'
sys.path.insert(0,str(AUDIT))
from align_structures import structure
A=pd.read_csv(AUDIT/'apistox_structures.csv').fillna('')
INITIAL=['cyclobutrifluram','isocycloseram','dimpropyridaz','acynonapyr','oxazosulfyl','fluhexafon','inpyrfluxam','spiropidion']
ADDITIONAL=['flupyrimin','tyclopyrazoflor','fluchlordiniliprole','benzpyrimoxan','flupentiofenox','fluxametamide','broflanilide','cyclopyrimorate']
def identity(name):
    u='https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+urllib.parse.quote(name)+'/property/InChIKey,InChI,IsomericSMILES/JSON'
    p=OUT/'cache'/(name+'_pubchem.json');p.parent.mkdir(exist_ok=True)
    old=OLD/(name+'_pubchem.json')
    try:
        if p.exists():j=json.loads(p.read_text(encoding='utf-8'))
        elif old.exists():j=json.loads(old.read_text(encoding='utf-8'))
        else:
            with urllib.request.urlopen(u,timeout=25) as r:j={'query_url':u,'response':json.loads(r.read())}
    except Exception as e:j={'query_url':u,'error':str(e)}
    p.write_text(json.dumps(j,indent=2),encoding='utf-8')
    rows=[]
    for v in j.get('response',j).get('PropertyTable',{}).get('Properties',[]):
        z=structure(v.get('SMILES',v.get('IsomericSMILES','')),v.get('InChI',''),v.get('InChIKey',''))
        x={'candidate':name,'candidate_origin':'previous8' if name in INITIAL else 'new8','PubChem_CID':v.get('CID'),'identity_url':u,**z}
        for k in ['inchikey','connectivity','parent_connectivity']:
            m=A[A[k].eq(z[k])] if z[k] else A.iloc[:0]
            x['overlap_'+k]=bool(len(m));x['matched_names_'+k]='|'.join(m.name)
        x['model_representation_note']='Chirality-insensitive ECFP permits connectivity representation; this does not prove exact material stereoisomer composition.'
        rows.append(x)
    return rows or [{'candidate':name,'structure_status':'unresolved','error':j.get('error','No properties')}]
SOURCES={
 'isocycloseram_pmra2025.pdf':'https://publications.gc.ca/collections/collection_2025/sc-hc/h113-9/H113-9-2025-11-eng.pdf',
 'dimpropyridaz_apvma.pdf':'https://www.apvma.gov.au/sites/default/files/publication/106036-public_release_summary_on_the_evaluation_of_the_new_active_constituent_dimpropyridaz_in_the_product_efficon_insecticide.pdf',
 'flupyrimin_frontiers2022.pdf':'https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2022.1019573/pdf',
 'flupyrimin_jafc2017.pdf':'https://irac-online.org/content/uploads/JAFC2017.pdf',
 'benzpyrimoxan_technical_guide.pdf':'https://www.nichino.co.jp/contents/000014547.pdf',
}
def source(item):
    name,u=item;p=OUT/'sources'/name;p.parent.mkdir(exist_ok=True);row={'file':str(p),'url':u}
    try:
        if not p.exists():
            with urllib.request.urlopen(u,timeout=30) as r:b=r.read()
            if not b.startswith(b'%PDF'):raise ValueError('Not PDF or empty response')
            p.write_bytes(b)
        row.update(status='downloaded',sha256=hashlib.sha256(p.read_bytes()).hexdigest())
        from pypdf import PdfReader
        pdf=PdfReader(str(p));p.with_suffix('.txt').write_text('\n'.join(f'\n=== PDF PAGE {i+1} ===\n'+(q.extract_text() or '') for i,q in enumerate(pdf.pages)),encoding='utf-8');row['pages']=len(pdf.pages)
    except Exception as e:row.update(status='unavailable',error=str(e))
    return row
if __name__=='__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:r=[z for rs in pool.map(identity,INITIAL+ADDITIONAL) for z in rs]
    F=pd.DataFrame(r);F.to_csv(OUT/'all_candidate_structures.csv',index=False,encoding='utf-8-sig')
    print(F[['candidate','PubChem_CID','overlap_parent_connectivity','matched_names_parent_connectivity']].to_string(index=False),flush=True)
    F[F.candidate_origin.eq('previous8')].to_csv(OUT/'initial8_structure_inference_panel.csv',index=False,encoding='utf-8-sig')
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:s=list(pool.map(source,SOURCES.items()))
    pd.DataFrame(s).to_csv(OUT/'source_manifest.csv',index=False,encoding='utf-8-sig')
    print(pd.DataFrame(s)[['file','status']].to_string(index=False))
