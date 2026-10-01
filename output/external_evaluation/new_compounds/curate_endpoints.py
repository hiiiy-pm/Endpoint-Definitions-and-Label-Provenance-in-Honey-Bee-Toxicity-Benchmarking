"""Document evidence adjudication made before any model prediction.

One row per extracted endpoint; never treat censored limits as exact LD50.
"""
from pathlib import Path
import json, urllib.request, urllib.parse, hashlib
import pandas as pd
OUT=Path(__file__).resolve().parent
ID=pd.read_csv(OUT/'candidate_structure_exclusion.csv').fillna('')
U={
 'acynonapyr':'https://publications.gc.ca/collections/collection_2026/sc-hc/h113-9/H113-9-2026-3-eng.pdf',
 'spiropidion':'https://publications.gc.ca/collections/collection_2026/sc-hc/h113-9/H113-9-2026-9-eng.pdf',
 'cyclobutrifluram':'https://publications.gc.ca/collections/collection_2025/sc-hc/h113-9/H113-9-2025-6-eng.pdf',
 'inpyrfluxam':'https://www.sumitomo-chem.co.jp/english/rd/report/2020E_1.pdf',
 'isocycloseram':'https://extranet.who.int/prequal/sites/default/files/doc_parts/WHOVC-SP_Isocycloseram_2025.1.pdf',
 'oxazosulfyl':'https://www.sumitomo-chem.co.jp/english/rd/report/2021E_2.pdf',
 'dimpropyridaz':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11981975/',
 'fluhexafon':''}
records=[]
def add(name,route,h,op,value,study,source_location,status,material,notes,unit='ug ai/bee'):
    row=dict(candidate=name,route=route,observation_h=h,qualifier=op,LD50_value=value,original_unit=unit,
             species='Apis mellifera' if name!='dimpropyridaz' else 'Honeybee (species not explicit in table)',
             life_stage='adult' if h else '',endpoint_name='LD50',study_id=study,source_url=U[name],
             source_location=source_location,evidence_status=status,test_material=material,notes=notes)
    if name=='dimpropyridaz':row['endpoint_name']='LC50 as printed; dose units require verification'
    for t in (1,11,100):
        row[f'label_le_{t}']=('0' if value>=t else 'unresolved') if op=='>' else ('1' if value<=t else '0') if op=='=' else 'unresolved'
    records.append(row)

add('acynonapyr','contact',48,'>',100,'PMRA 3328790 (2016)','Table 24 printed pp106-107/PDF112-113; reference printed p149/PDF155',
    'regulatory_fulltext_verified','technical active ingredient; bee-study purity 99%',
    'Full regulator PDF text reviewed via web. Local host returned empty HTTP202; source download unavailable. Study ID confirmed to cover oral AND contact. CAS1332838-17-1. No precise LD50. Regulator specifies3-endo; PubChem structure has unspecified attached-center stereochemistry. CAS agreement does not establish complete material-level stereochemical identity; connectivity exclusion remains valid.')
add('acynonapyr','oral',48,'>',98.8,'PMRA 3328790 (2016)','Table24 printed p106/PDF112','regulatory_fulltext_verified','technical active ingredient; bee-study purity 99%',
    'Cannot label at 100 threshold; limit98.8 is below100. Do not round to100. Regulator specifies3-endo; PubChem stereochemical completeness needs verification; CAS match alone is insufficient for full material identity.')
add('spiropidion','contact',48,'>',200,'PMRA 3436778 (2018)','Table19 printed p92/PDF98; reference printed p153/PDF159',
    'regulatory_fulltext_verified','active ingredient SYN546330; bee-study batch purity not shown',
    'Parent active ingredient explicitly distinguished from formulation A20262B and metabolite SYN547305. General technical purity98.3% must NOT be treated as study batch purity. CAS1229023-00-0.')
add('spiropidion','oral',48,'>',100,'PMRA 3436778 (2018)','Table19 printed p92/PDF98',
    'regulatory_fulltext_verified','active ingredient SYN546330; bee-study batch purity not shown','Right-censored; not an exact value of100.')
add('inpyrfluxam','contact',48,'>',100,'Not reported in developer report; EPA index points to MRID49706023','Table9 printed/PDF p12; technical material clarified p13',
    'developer_fulltext_verified','technical product; batch purity not reported',
    'Developer primary report downloaded and table rendered/read. EPA precise study match remains index-only, not counted as second independent experiment. Unit printed ug/bee; active technical material stated. Adult stage is not explicitly printed; adult-test inference requires underlying study verification. No original studyID or batch purity in this report.',unit='ug/bee')
add('inpyrfluxam','oral',48,'>',111.3,'Not reported in developer report','Table9 printed/PDF p12; technical material clarified p13',
    'developer_fulltext_verified','technical product; batch purity not reported','APVMA chemistry source CAS1352994-67-2; not a formulation. Adult stage is not explicitly printed; adult-test inference requires underlying study verification. No original studyID or batch purity in this report.',unit='ug/bee')
add('cyclobutrifluram','contact',48,'>',200,'PMRA3273332 (2020); EPA MRID51459447 index-only','Table22 printed p126/PDF131; identity printed pp9-10/PDF14-15',
    'endpoint_verified_identity_hold','enantiomer mixture; registered technical purity85% is not bee-batch purity',
    'Registered CAS1460292-16-3 mixture80-100%(1S,2S) and0-20%(1R,2R); previous PubChem name lookup is one stereoisomer. Connectivity excludes training overlap but full InChIKey not an identity proof for tested mixture. Hold strict inclusion.')
add('cyclobutrifluram','oral',48,'>',72.2,'PMRA3273332 (2020)','Table22 printed p126/PDF131','endpoint_verified_identity_hold','enantiomer mixture',
    'PMRA reports72.2; EPA index reports72.23. Keep primary-read72.2; do not count as independent studies. 100 threshold unresolved.')
add('isocycloseram','contact',96,'=',0.26,'VV-466340; Report S16-02527 (Kling2016)','WHO2025 Table6 printed/PDF p34; reference p41','verified_nonmatching_time','technical active; purity96.9%',
    '96h cannot be relabelled48h. WHO48h0.35/>10 values are Bombus terrestris, not Apis. MN0.072/0.18 values refer to a degradate, not parent.',unit='ug/bee')
add('isocycloseram','oral',72,'=',0.28,'VV-466340; Report S16-02527 (Kling2016)','WHO2025 Table6 printed/PDF p34; reference p41','verified_nonmatching_time','technical active; purity96.9%','72h cannot be relabelled48h.',unit='ug/bee')
add('oxazosulfyl','contact',72,'=',0.077,'Not reported','Developer2021 Table8 printed/PDF p13','verified_nonmatching_time','technical grade active ingredient','Downloaded original developer PDF;72h not48h.',unit='ug/bee')
add('oxazosulfyl','oral',96,'=',0.015,'Not reported','Developer2021 Table8 printed/PDF p13','verified_nonmatching_time','technical grade active ingredient','96h not48h.',unit='ug/bee')
add('dimpropyridaz','contact',48,'>',50.3,'Not traced','Pest Management Science2025 Table8','index_only_endpoint_name_unresolved','not established from full study',
    'Indexed table calls endpoint LC50 although units are per-bee dose. PMC full text inaccessible at collection time; species/study-material/ID require regulator report. Threshold100 unresolved.')
add('dimpropyridaz','oral',48,'>',43.3,'Not traced','Pest Management Science2025 Table8','index_only_endpoint_name_unresolved','not established from full study','Same source limitations; threshold100 unresolved.')
add('fluhexafon','',None,'',float('nan'),'Not found','','no_primary_endpoint_found','','No suitable primary bee endpoint found within bounded search. Retained in candidate ledger, never excluded by model score.')
E=pd.DataFrame(records).merge(ID,on='candidate',validate='many_to_one')
E.loc[E.candidate.eq('inpyrfluxam'),'life_stage']='not_explicit_in_source'
E['complete_material_identity_confirmed']=False
E['candidate_48h_panel']=E.evidence_status.isin(['regulatory_fulltext_verified','developer_fulltext_verified']) & E.observation_h.eq(48)
E.to_csv(OUT/'candidate_endpoints.csv',index=False,encoding='utf-8-sig')
E[E.candidate_48h_panel].to_csv(OUT/'reviewed_48h_external_candidate_panel.csv',index=False,encoding='utf-8-sig')
E[E.evidence_status.eq('verified_nonmatching_time')].to_csv(OUT/'extended_time_candidates_not_48h.csv',index=False,encoding='utf-8-sig')
coverage=[]
for route in ['contact','oral']:
    panel=E[E.candidate_48h_panel&E.route.eq(route)]
    for t in (1,11,100):
        counts=panel[f'label_le_{t}'].value_counts()
        coverage.append(dict(route=route,threshold=t,n_compounds=len(panel),positive=int(counts.get('1',0)),negative=int(counts.get('0',0)),unresolved=int(counts.get('unresolved',0))))
pd.DataFrame(coverage).to_csv(OUT/'panel_class_coverage.csv',index=False,encoding='utf-8-sig')
pd.DataFrame([dict(panel='48h_contact',le1=0,gt1_le11=0,gt11_le100=0,gt100=3,unresolved=0),dict(panel='48h_oral',le1=0,gt1_le11=0,gt11_le100=0,gt100=2,unresolved=1)]).to_csv(OUT/'panel_four_tier_coverage.csv',index=False,encoding='utf-8-sig')
checks=[]
def check(name,result,detail):checks.append(dict(check=name,passed=bool(result),detail=detail))
check('fixed_candidate_ledger',set(E.candidate)==set(ID.candidate),'8 candidates retained including no-data and rejected candidates.')
check('structure_keys_complete',ID.inchikey.str.len().eq(27).all(),'RDKit identities from cached PubChem; fullkey,connectivity,parent checked.')
check('training_exclusion_all_three',not ID[['overlap_inchikey','overlap_connectivity','overlap_parent_connectivity']].any().any(),'All8 connectivity structures absent from frozen1035; mixture caveat retained.')
check('isomers_not_claimed_equivalent',not E.loc[E.candidate.eq('cyclobutrifluram'),'candidate_48h_panel'].any(),'Cyclobutrifluram enantiomer mixture held despite no connectivity overlap.')
check('time_not_coerced',not E.loc[E.candidate.isin(['isocycloseram','oxazosulfyl']),'candidate_48h_panel'].any(),'72/96h remain72/96h.')
check('limits_not_exact',(E.loc[E.candidate_48h_panel,'qualifier']=='>').all(),'Every eligible panel endpoint is right-censored; exact regression validation unavailable.')
check('98.8_at100_unresolved',E.loc[E.candidate.eq('acynonapyr')&E.route.eq('oral'),'label_le_100'].iloc[0]=='unresolved','No rounding across decision threshold.')
check('one_class_warning',(E.loc[E.candidate_48h_panel&E.route.eq('contact'),'label_le_100']=='0').all(),'3 contact chemicals all negative at1/11/100; AUROC undefined, not a sufficient standalone set.')
check('WHO_original_unit_preserved',E.loc[E.candidate.eq('isocycloseram'),'original_unit'].eq('ug/bee').all(),'WHO printsug/bee; do not insert ai into original unit.')
check('inpyrfluxam_stage_not_invented',E.loc[E.candidate.eq('inpyrfluxam'),'life_stage'].eq('not_explicit_in_source').all(),'48h candidate retained; adult stage needs original study confirmation.')
pd.DataFrame(checks).to_csv(OUT/'first_round_checks.csv',index=False,encoding='utf-8-sig')
assert all(c['passed'] for c in checks)
print(E[['candidate','route','observation_h','qualifier','LD50_value','candidate_48h_panel']].to_string(index=False))
print('First-round automated checks:',len(checks),'passed')
