"""Check regulator CAS against PubChem structure independently of name lookup."""
from pathlib import Path
import pandas as pd
from collect_candidates import identity
OUT=Path(__file__).resolve().parent
N=pd.read_csv(OUT/'candidate_structure_exclusion.csv').set_index('candidate')
rows=[]
for name,cas in [('acynonapyr','1332838-17-1'),('spiropidion','1229023-00-0'),('inpyrfluxam','1352994-67-2'),('cyclobutrifluram','1460292-16-3')]:
    for r in identity(cas):
        r['candidate']=name;r['regulator_CAS']=cas
        r['name_CAS_fullkey_agree']=r.get('inchikey','')==N.loc[name,'inchikey']
        r['name_CAS_connectivity_agree']=r.get('connectivity','')==N.loc[name,'connectivity']
        r['cas_is_overlap_key']=False
        r['mixture_caution']='CAS query alone does not prove stereoisomer-ratio representation' if name=='cyclobutrifluram' else ''
        rows.append(r)
F=pd.DataFrame(rows)
F.to_csv(OUT/'regulatory_cas_crosscheck.csv',index=False,encoding='utf-8-sig')
print(F[['candidate','regulator_CAS','PubChem_CID','name_CAS_fullkey_agree','name_CAS_connectivity_agree']].to_string(index=False))
