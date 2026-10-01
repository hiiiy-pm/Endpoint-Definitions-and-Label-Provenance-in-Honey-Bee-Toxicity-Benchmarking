"""Check scientific boundary logic, source cells, counts and source integrity."""
from pathlib import Path
import hashlib
import json
import platform
import sys
import importlib.metadata

import openpyxl
import pandas as pd
from audit_external_data import label

OUT=Path(__file__).resolve().parent
BASE=(Path(__file__).resolve().parents[2] / 'data/external_raw')
checks={}
cases=[(100,'>',None,'',100,0),(100,'>=',None,'',100,-1),
       (100,'',None,'',100,1),(11,'>',None,'',100,-1),
       (11,'>',None,'',11,0),(1,'<',None,'',1,1),
       (None,'',1,'<=',1,1),(0.1,'>',10,'<',11,1),
       (1,'ca.',None,'',11,-1)]
for lo,lq,hi,hq,t,expected in cases:
    assert label(lo,lq,hi,hq,t)==expected,(lo,lq,hi,hq,t)
checks['interval_boundary_cases']=len(cases)
inv=json.loads((BASE/'INPUT_MANIFEST.json').read_text(encoding='utf-8'))
for item in inv:
    assert hashlib.sha256((BASE/item['file']).read_bytes()).hexdigest()==item['sha256']
checks['source_hashes_verified']=len(inv)
o=pd.read_csv(OUT/'oft_acute_endpoints_aligned.csv')
p=pd.read_csv(OUT/'plos_adult_acute_endpoints_aligned.csv')
assert o.result_uuid.is_unique
assert len(o)==751 and o.ld50_screen.sum()==738
assert len(p)==319 and p.substance_id.nunique()==173
checks['source_and_endpoint_count_invariants']=True
summary=pd.read_csv(OUT/'overlap_summary.csv')
for key in ['inchikey','connectivity','parent_connectivity']:
    assert (summary['n_'+key]==summary['overlap_'+key]+summary['unmatched_'+key]).all()
checks['overlap_denominators_reconcile']=True
for table in [o,p]:
    known=table[[f'y_{t}' for t in [1,11,100]]].ge(0).all(axis=1)
    assert (table.loc[known,'y_1']<=table.loc[known,'y_11']).all()
    assert (table.loc[known,'y_11']<=table.loc[known,'y_100']).all()
    assert not (table.conservative_unseen_candidate & table.inchikey.isna()).any()
checks['nested_labels_and_unresolved_identity_policy']=True
wb=openpyxl.load_workbook(BASE/'journal.pone.0265962.s001.xlsx',read_only=True,data_only=True)
ws=wb.active
for cell in ['F8','F9','F11','G11','F15','F201']:
    row=p[p.source_cell.eq(cell)].iloc[0]
    assert str(ws[cell].value).strip()==str(row.raw_value).strip(),cell
wb.close()
checks['S1_source_cells_verified']=6
assert str(p[(p['name']=='Abamectin') & (p.route=='Contact')].iloc[0].MRID)=='49226803'
assert str(p[(p['name']=='Acetamiprid') & (p.route=='Oral')].iloc[0].MRID)=='44651874'
checks['S4_spot_checks']=2
raw=pd.read_pickle(OUT/'cache/END_STUDY_REC.TerrestEcotox.pkl').set_index('excel_row')
for _,row in o.sample(12,random_state=203).iterrows():
    source=raw.loc[int(str(row.source_excel_rows).split('|')[0])]
    assert source['ResultsAndDiscussion.EffectConcentrations.UUID']==row.result_uuid
    assert source['Parent UUID']==row.substance_id
checks['OFT_source_rows_verified']=12
u=pd.read_csv(OUT/'excluded_all_parsed.csv')
assert len(u)==1389 and u.source_block.nunique()==23
g=pd.read_csv(OUT/'excluded_uncertain_route_groups.csv')
assert g['size'].sum()==745 and g.CAS.nunique()==72 and g.upstream_unspecified_rule.sum()==79
checks['excluded_blocks_and_groups_reconcile']=True
con=pd.read_csv(OUT/'label_concordance_summary.csv')
assert (con.n==con.agreements+con.disagreements).all()
assert con[con.cohort.eq('all_three_determinate')].groupby('source').n.nunique().eq(1).all()
checks['fixed_cohort_disagreement_denominators_reconcile']=True
checks['environment']={'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
    **{k:importlib.metadata.version(k) for k in ['pandas','numpy','openpyxl','rdkit']}}
(OUT/'verification.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(checks,ensure_ascii=False,indent=2))
