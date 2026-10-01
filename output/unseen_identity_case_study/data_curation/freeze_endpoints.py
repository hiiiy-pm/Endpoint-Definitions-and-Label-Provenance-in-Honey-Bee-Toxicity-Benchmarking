"""Freeze evidence-based endpoints without reading any model prediction."""
from pathlib import Path
import pandas as pd,json,hashlib
OUT=Path(__file__).resolve().parent
OLD=OUT.parents[1]/'external_evaluation/new_compounds'
S=pd.read_csv(OUT/'all_candidate_structures.csv').fillna('')
E=pd.read_csv(OLD/'candidate_endpoints.csv').fillna('')
cols=['candidate','route','observation_h','qualifier','LD50_value','original_unit','species','life_stage','endpoint_name','study_id','source_url','source_location','evidence_status','test_material','notes']
E=E[cols].copy()
E['new_primary_48h_eligible']=False
E['eligibility_reason']='Retained previous record; not eligible unless explicitly updated.'
for name in ['acynonapyr','spiropidion','cyclobutrifluram']:
    mask=E.candidate.eq(name)
    E.loc[mask,'new_primary_48h_eligible']=True
    E.loc[mask,'eligibility_reason']='Regulatory dose LD50;48h;Apis;adult/nonlarval context;single connectivity. Exact material stereochemistry remains unconfirmed.'
E.loc[E.candidate.eq('cyclobutrifluram'),'evidence_status']='regulatory_fulltext_verified_connectivity_unit'
E.loc[E.candidate.eq('cyclobutrifluram'),'notes']+=' Same-connectivity enantiomer mixture accepted only for exploratory chirality-insensitive model; full InChIKey not exact material identity. Adult interpreted from regulator separate larval/brood rows and adult-risk endpoint use.'
E.loc[E.candidate.eq('cyclobutrifluram'),'notes']=E.loc[E.candidate.eq('cyclobutrifluram'),'notes'].str.replace('Hold strict inclusion.','Previous20260922identityhold is replaced by20260923 pre-prediction connectivity-unit policy.',regex=False)
E.loc[E.candidate.eq('spiropidion'),'notes']+=' Adult-stage interpretation from regulatory table separating acute honeybee rows from larval rows; study-specific batch purity absent.'
E.loc[E.candidate.eq('inpyrfluxam'),'eligibility_reason']='48h and dose verified, but life stage not explicit in developer report; withheld from primary.'
E.loc[E.candidate.isin(['isocycloseram','oxazosulfyl']),'eligibility_reason']='Keep72/96h as extended-time only; never relabel48h.'
E.loc[E.candidate.eq('dimpropyridaz'),'eligibility_reason']='APVMA full text supports adult-bee dose LD50 bounds but combines2 species and omitsobservation time; developer table saysLC50. No48h primary assignment.'
E.loc[E.candidate.eq('dimpropyridaz'),'notes']+=' APVMA report p32/PDF36 gives lowest adult contactLD50>50 and oral>42 across2 species; it does not independently establish species-specific48h50.3/43.3 values.'
rows=[]
def add(name,route,h,q,val,unit,study,url,loc,status,material,notes,eligible,stage='adult',species='Apis mellifera',endpoint='LD50'):
    rows.append(dict(candidate=name,route=route,observation_h=h,qualifier=q,LD50_value=val,original_unit=unit,study_id=study,source_url=url,source_location=loc,evidence_status=status,test_material=material,notes=notes,new_primary_48h_eligible=eligible,eligibility_reason='Regulatory primary dose endpoint verified' if eligible else notes,life_stage=stage,species=species,endpoint_name=endpoint))
pmra='https://publications.gc.ca/collections/collection_2025/sc-hc/h113-9/H113-9-2025-11-eng.pdf'
add('isocycloseram','oral',48,'=',0.29,'ug ai/bee','PMRA3246199;2016;SYN547407',pmra,'Table17 printedp133/PDF138; source lines7283-7294;study reference printedp180/PDF185','regulatory_fulltext_verified','Isocycloseram Technical; table grouped purity96.9-98.8%w/w',
    'Endpoint column explicitly48h0.29 although study observation column says72h adult. Parent active,notSYN549106;Apis,notBombus. WHO72h0.28 likely different regulatory rounding/interpretation of2016 study; not independent replication.',True)
add('isocycloseram','oral',72,'=',0.29,'ug ai/bee','PMRA3246199;2016;SYN547407',pmra,'Table17 printedp133/PDF138','regulatory_fulltext_verified_extended','Isocycloseram Technical; grouped purity96.9-98.8%',
    'Extended only; same likely2016 study asWHO72h0.28. Do not double count as independent observation.',False)
flux='https://www.maff.go.jp/j/nouyaku/n_sinsa/attach/pdf/index-38.pdf'
for route,val in [('contact',2.68),('oral',0.792)]:
    add('fluxametamide',route,48,'=',val,'ug/bee (source: µg/頭)','NC-515 Acute Toxicity to Honey Bees;2014;HuntingdonLifeSciencesONE0224(E0701)',flux,'Table2.6-11 printedp106/PDF113; bibliography II.2.6.3.1','regulatory_fulltext_and_visual_verified','Technical active 原体;batch purity not printed inbee table',
        'Apis mellifera成虫 explicit;10bees/group3replicates. Direct dose experiment,noLC50 conversion. Table2.6-12 below is separate10%EC field study and not used. Unit does not printai; noaddedai qualifier. MaterialCAS928783-29-3;racemic5RS;connectivity-level model unit.',True)
flup='https://www.maff.go.jp/j/council/sizai/nouyaku/attach/pdf/mitubati_17-7.pdf'
for route,val,num in [('contact',100,'13912.6119;KCI200136'),('oral',500,'13912.6118;KCI200140')]:
    add('flupentiofenox',route,48,'>',val,'ug ai/bee','SmithersUSA;2020;Report'+num,flup,'Tables2-3 printedpp8-9/PDF11-12;references PDF23','regulatory_fulltext_and_visual_verified','Technical active 原体;bee-batch purity not printed',
        'AdultApis explicit;OECD214contact/213oral;10bees/group3replicates. Retaingreater-than from actualstudy tables; regulatoryrisk metric later dropsqualifier butisnotexactLD50. CAS1472050-04-6.',True)
guide='https://www.nichino.co.jp/contents/000014547.pdf'
for route in ['contact','oral']:
    add('benzpyrimoxan',route,48,'>',100,'ug ai/bee','Notgiven',guide,'Technicalguide PDF6 bottomtable','developer_visual_verified_material_hold','BPX;testmaterial technicalversusformulation notexplicit',
        'AdultApis48h explicitanddoseunit clear;combinedLD50/LC50 tableheading. MaterialandstudyID unresolved;heldoutsideprimary.',False)
add('fluxametamide','oral','', '=',1.083,'ug/adult','Song2024;DOI10.1016/j.pestbp.2024.106109','https://pubmed.ncbi.nlm.nih.gov/39277414/','Abstract','abstract_only_time_unverified','Notverified',
    'Potentialdifferent-study estimate;observationtimeandspecies/batchnotverified. Do notcherry-pickagainstregulatory0.792.',False)
add('fluxametamide','contact','', '=',0.25,'ug ai/bee','Liu2026;DOI10.1002/ps.70708','https://pubmed.ncbi.nlm.nih.gov/41793051/','Abstract','abstract_only_time_unverified','Notverified',
    'Potentialdifferent-study estimate;observationtimeunverified. NotprimaryandnotaveragedwithMAFF2.68.',False)
for name in ['tyclopyrazoflor','fluchlordiniliprole','cyclopyrimorate']:
    add(name,'','','',None,'','Notfound','','boundedsearch','no_primary_endpoint_found','Unknown','No verified primaryApis doseLD50 withinboundedsearch.',False,stage='',species='')
for name in ['flupyrimin','broflanilide']:
    add(name,'','','',None,'','Excludedbystructure','','all_candidate_structures.csv','excluded_training_overlap','','Connectivity/parent overlapsfrozen1035;notexternalnewcompound.',False,stage='',species='')
E=pd.concat([E,pd.DataFrame(rows)],ignore_index=True)
E['compound_id']=E.candidate
E['life_stage_evidence']='explicit_in_source' 
E.loc[E.candidate.isin(['spiropidion','cyclobutrifluram']),'life_stage_evidence']='adult_from_regulatory_context;not_explicit_in_endpoint_cell'
E.loc[E.candidate.eq('inpyrfluxam'),'life_stage_evidence']='not_explicit_in_source'
E['endpoint_id']=[f'EXT20260923_{i+1:03d}' for i in range(len(E))]
E=E.merge(S[['candidate','inchikey','connectivity','parent_connectivity','overlap_inchikey','overlap_connectivity','overlap_parent_connectivity']],on='candidate',validate='many_to_one')
E['new_primary_48h_eligible']=E.new_primary_48h_eligible & ~E.overlap_parent_connectivity
def label(r,t):
    try:v=float(r.LD50_value)
    except:return ''
    if r.qualifier=='=':return int(v<=t)
    if r.qualifier=='>' and v>=t:return 0
    return ''
for t in [1,11,100]:E[f'label_le_{t}']=[label(r,t) for r in E.itertuples()]
E.to_csv(OUT/'frozen_endpoint_ledger.csv',index=False,encoding='utf-8-sig')
P=E[E.new_primary_48h_eligible].copy()
P.to_csv(OUT/'frozen_scoreable_48h_endpoints.csv',index=False,encoding='utf-8-sig')
P[P.route.eq('contact')].to_csv(OUT/'frozen_primary_48h_contact.csv',index=False,encoding='utf-8-sig')
P[P.route.eq('oral')].to_csv(OUT/'frozen_secondary_48h_oral.csv',index=False,encoding='utf-8-sig')
P[P.life_stage_evidence.eq('explicit_in_source')].to_csv(OUT/'conservative_explicit_adult_48h_subset.csv',index=False,encoding='utf-8-sig')
E[E.observation_h.isin([24,72,96,24.0,72.0,96.0])].to_csv(OUT/'extended_times_not_primary.csv',index=False,encoding='utf-8-sig')
cov=[]
for route,G in P.groupby('route'):
    d=dict(route=route,n_compounds=G.candidate.nunique(),tier_le1=0,tier_1_11=0,tier_11_100=0,tier_gt100=0,tier_unresolved=0)
    for r in G.itertuples():
        v=float(r.LD50_value)
        k='tier_le1' if r.qualifier=='=' and v<=1 else 'tier_1_11' if r.qualifier=='=' and v<=11 else 'tier_11_100' if r.qualifier=='=' and v<=100 else 'tier_gt100' if v>100 or r.qualifier=='>' and v==100 else 'tier_unresolved'
        d[k]+=1
    for t in [1,11,100]:
        d[f'positive_{t}']=sum(G[f'label_le_{t}'].eq(1));d[f'negative_{t}']=sum(G[f'label_le_{t}'].eq(0));d[f'unresolved_{t}']=sum(G[f'label_le_{t}'].eq(''))
    cov.append(d)
pd.DataFrame(cov).to_csv(OUT/'frozen_panel_coverage.csv',index=False,encoding='utf-8-sig')
added=S[S.candidate.isin(['fluxametamide','flupentiofenox'])].copy();added['compound_id']=added.candidate;added['name']=added.candidate;added['SMILES']=added.canonical_smiles
added[['compound_id','name','SMILES','inchikey','connectivity','parent_connectivity','PubChem_CID']].to_csv(OUT/'frozen_added2_inference_structures.csv',index=False,encoding='utf-8-sig')
manifest={'no_prediction_files_read':True,'candidate_count':16,'primary_endpoint_rows':len(P),'files':{}}
for name in ['frozen_endpoint_ledger.csv','frozen_scoreable_48h_endpoints.csv','frozen_primary_48h_contact.csv','frozen_secondary_48h_oral.csv','frozen_added2_inference_structures.csv']:
    manifest['files'][name]=hashlib.sha256((OUT/name).read_bytes()).hexdigest()
(OUT/'endpoint_freeze_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
checks={
 'all16retained':E.candidate.nunique()==16,
 'primary48h_only':P.observation_h.astype(float).eq(48).all(),
 'Apis_only':P.species.eq('Apis mellifera').all(),
 'none_in_training':not P.overlap_parent_connectivity.any(),
 'no_unknown_stage_primary':not P.life_stage.eq('not_explicit_in_source').any(),
 'positive_from_parent_notdegradate':len(P[P.candidate.eq('isocycloseram')&P.route.eq('oral')&P.LD50_value.eq(0.29)])==1,
 'flux_contact_middle':len(P[P.candidate.eq('fluxametamide')&P.route.eq('contact')&P.label_le_1.eq(0)&P.label_le_11.eq(1)&P.label_le_100.eq(1)])==1,
 'flupentio_limits_preserved':P.loc[P.candidate.eq('flupentiofenox'),'qualifier'].eq('>').all(),
 'inpyrfluxam_held':not P.candidate.eq('inpyrfluxam').any(),
 'dimpropyridaz_held':not P.candidate.eq('dimpropyridaz').any(),
 'unique_compound_route':not P.duplicated(['candidate','route']).any(),
}
pd.DataFrame([dict(check=k,passed=v) for k,v in checks.items()]).to_csv(OUT/'round1_data_checks.csv',index=False,encoding='utf-8-sig')
assert all(checks.values()),checks
print(pd.DataFrame(cov).to_string(index=False));print('checks',checks)
