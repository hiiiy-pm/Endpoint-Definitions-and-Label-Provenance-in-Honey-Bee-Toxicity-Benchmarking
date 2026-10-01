"""Reproducible feasibility audit, with no model fitting or source changes.

Dependencies: pandas, numpy, openpyxl, RDKit. Run inspect_inputs.py,
extract_tables.py, align_structures.py, fetch_provenance.py first.
"""
from pathlib import Path
import hashlib
import json
import re
from collections import Counter

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from align_structures import norm, s

OUT=Path(__file__).resolve().parent
CACHE=OUT/'cache'
RDLogger.DisableLog('rdApp.*')
THRESHOLDS=(1,11,100)
API=pd.read_csv(OUT/'apistox_structures.csv').fillna('')
KEYS=('inchikey','connectivity','parent_connectivity')
ASETS={k:set(API[k])-{''} for k in KEYS}
ACAS={norm(c) for c in API.CAS if s(c)}
ANAMES={norm(re.sub(r'\s*\(Ref:.*','',n,flags=re.I)) for n in API.name}

def save(df,name):
    df.to_csv(OUT/name,index=False,encoding='utf-8-sig')

def finite(x):
    try:
        v=float(x)
        return v if np.isfinite(v) else None
    except (ValueError,TypeError):
        return None

def label(lo,lq,hi,hq,t):
    """y=1 means LD50 <= t. -1 means not identifiable from the interval."""
    lo,hi=finite(lo),finite(hi)
    lq,hq=s(lq),s(hq)
    if lo is not None and hi is None and lq in ('','='):
        return int(lo<=t)
    if lo is not None and lq in ('<','<=') and lo<=t:
        return 1
    if hi is not None and hq in ('<','<=','') and hi<=t:
        return 1
    if lo is not None and lq in ('>','>=') and (lo>t or (lo==t and lq=='>')):
        return 0
    return -1

def qualify(df):
    for k in KEYS:
        df['overlap_'+k]=df[k].fillna('').isin(ASETS[k]) & df[k].fillna('').ne('')
    df['same_CAS_in_apistox']=df['CAS'].map(lambda x:norm(x) in ACAS if s(x) else False)
    df['same_name_in_apistox']=df['name'].map(lambda x:norm(x) in ANAMES if s(x) else False)
    df['identity_conflict_or_same_chemical']=df['same_CAS_in_apistox']|df['same_name_in_apistox']
    df['conservative_unseen_candidate']=(df['inchikey'].fillna('').ne('') &
        df['parent_connectivity'].fillna('').ne('') & ~df['overlap_inchikey'] &
        ~df['overlap_connectivity'] & ~df['overlap_parent_connectivity'] &
        ~df['identity_conflict_or_same_chemical'] & ~df['inchi_smiles_disagree'])
    return df

def counts(d,source,subset,route):
    ans={'source':source,'subset':subset,'route':route,'records':len(d),
         'source_substances':d['substance_id'].nunique(),
         'unresolved_substances':d.loc[d.inchikey.fillna('').eq(''),'substance_id'].nunique()}
    for key in KEYS:
        vals=set(d[key].dropna())-{''}
        ans['n_'+key]=len(vals)
        ans['overlap_'+key]=len(vals&ASETS[key])
        ans['unmatched_'+key]=len(vals-ASETS[key])
    ans['conservative_unseen_full_keys']=d.loc[d.conservative_unseen_candidate,'inchikey'].nunique()
    ans['exact_records']=int(d['value_kind'].eq('exact').sum())
    return ans

def threshold_tables(d,source,subset):
    summary=[]
    for route,r in d.groupby('route'):
        for t in THRESHOLDS:
            labels=[]
            for key,g in r[r.inchikey.fillna('').ne('')].groupby('inchikey'):
                vals=set(g[f'y_{t}'])
                status='conflict' if 0 in vals and 1 in vals else 'unknown' if -1 in vals else 'positive' if vals=={1} else 'negative'
                labels.append(status)
            summary.append({'source':source,'subset':subset,'route':route,'threshold_ug_bee':t,
                'records':len(r),'record_positive':int(r[f'y_{t}'].eq(1).sum()),
                'record_negative':int(r[f'y_{t}'].eq(0).sum()),'record_unknown':int(r[f'y_{t}'].eq(-1).sum()),
                'unique_structures':len(labels),**{k:labels.count(k) for k in ['positive','negative','unknown','conflict']}})
    return summary

def oft():
    raw=pd.read_pickle(CACHE/'END_STUDY_REC.TerrestEcotox.pkl')
    bee=raw[raw['MaterialsAndMethods.TestOrganisms.TestOrganismsSpecies'].eq('Apis mellifera')].copy()
    save(bee,'oft_bee_source_rows.csv')
    selected=bee[bee['AdministrativeData.Endpoint'].isin(['toxicity to bees: acute contact','toxicity to bees: acute oral'])].copy()
    pre='ResultsAndDiscussion.EffectConcentrations.'
    duplicate_check_cols=[c for c in selected if c.startswith(pre)]
    duplicate_variants=selected.groupby(pre+'UUID')[duplicate_check_cols].nunique(dropna=False).max().max()
    assert duplicate_variants==1, 'Duplicated result UUID has conflicting result fields'
    excel_map=selected.groupby(pre+'UUID').excel_row.apply(lambda x:'|'.join(map(str,x)))
    ref_map=selected.groupby(pre+'UUID')['DataSource.Reference'].apply(lambda x:'|'.join(sorted(set(x.dropna()))))
    d=selected.drop_duplicates(pre+'UUID').copy()
    subs=pd.read_pickle(CACHE/'SUB.pkl').set_index('Document UUID')
    refs=pd.read_pickle(CACHE/'oft_reference_structures.pkl').set_index('Document UUID')
    lit=pd.read_pickle(CACHE/'LIT.pkl').set_index('Document UUID')
    z=pd.DataFrame(index=d.index)
    z['result_uuid']=d[pre+'UUID']
    z['document_uuid']=d['Document UUID']
    z['source_excel_rows']=z.result_uuid.map(excel_map)
    z['substance_id']=d['Parent UUID']
    z['name']=z.substance_id.map(subs.ChemicalName)
    z['composition']=z.substance_id.map(subs['TypeOfSubstance.Composition'])
    z['reference_uuid']=z.substance_id.map(subs['ReferenceSubstance.ReferenceSubstance'])
    for key in ['inchikey','connectivity','parent_connectivity','canonical_smiles','inchi_smiles_disagree','fragments']:
        z[key]=z.reference_uuid.map(refs[key])
    z['CAS']=z.reference_uuid.map(refs['Inventory.CASNumber'])
    z['reference_smiles']=z.reference_uuid.map(refs['MolecularStructuralInfo.SmilesNotation'])
    z['reference_inchi']=z.reference_uuid.map(refs['MolecularStructuralInfo.InChl'])
    z['route']=d['AdministrativeData.Endpoint'].map({'toxicity to bees: acute contact':'Contact','toxicity to bees: acute oral':'Oral'})
    z['endpoint']=d[pre+'Endpoint']
    z['endpoint_other']=d[pre+'Endpoint.Other']
    z['unit_raw']=d[pre+'EffectConc.Unit']
    z['unit_other_raw']=d[pre+'EffectConc.Unit.Other']
    z['unit']=z.unit_raw.where(z.unit_raw.ne('other:'),z.unit_other_raw)
    z['unit_mapping']=z.unit.map(lambda v:'direct_ug_bee' if s(v) in ('µg/bee','μg/bee') else 'species_based_piece_to_bee' if s(v) in ('µg/piece','μg/piece') else 'mg_bee_to_ug_bee' if s(v)=='mg/bee' else 'ineligible_unit')
    factor=z.unit_mapping.map({'direct_ug_bee':1.,'species_based_piece_to_bee':1.,'mg_bee_to_ug_bee':1000.})
    z['lower_ug_bee']=pd.to_numeric(d[pre+'EffectConc.lowerValue'],errors='coerce')*factor
    z['lower_qualifier']=d[pre+'EffectConc.lowerQualifier'].fillna('')
    z['upper_ug_bee']=pd.to_numeric(d[pre+'EffectConc.upperValue'],errors='coerce')*factor
    z['upper_qualifier']=d[pre+'EffectConc.upperQualifier'].fillna('')
    z['value_kind']=np.where(z.lower_ug_bee.notna() & z.lower_qualifier.isin(['','=']) & z.upper_ug_bee.isna(),'exact',np.where(z.lower_ug_bee.notna()|z.upper_ug_bee.notna(),'censored_or_approximate','missing_or_ineligible'))
    z['duration_value']=d[pre+'Duration.Value'].combine_first(d['MaterialsAndMethods.StudyDesign.TotalExposureDuration.Value'])
    z['duration_unit']=d[pre+'Duration.Unit'].combine_first(d['MaterialsAndMethods.StudyDesign.TotalExposureDuration.Unit'])
    z['duration_hours']=pd.to_numeric(z.duration_value,errors='coerce')*z.duration_unit.map({'h':1,'d':24})
    z['test_material']=d['MaterialsAndMethods.TestMaterials.SpecificDetailsOnTestMaterialUsedForTheStudy']
    z['dose_basis']=d[pre+'ConcBasedOn']
    z['dose_basis_other']=d[pre+'ConcBasedOn.Other']
    z['basis_for_effect']=d[pre+'BasisForEffect']
    z['material_review_needed']=z.test_material.notna()
    z['literature_uuids']=z.result_uuid.map(ref_map)
    z['source_dois']=z.literature_uuids.map(lambda v:'|'.join(sorted(set(s(lit.loc[k,'GeneralInfo.Source']) for k in s(v).split('|') if k in lit.index))))
    z['source_years']=z.literature_uuids.map(lambda v:'|'.join(sorted(set(s(lit.loc[k,'GeneralInfo.ReferenceYear']) for k in s(v).split('|') if k in lit.index))))
    z=qualify(z)
    for t in THRESHOLDS:
        z[f'y_{t}']=[label(r.lower_ug_bee,r.lower_qualifier,r.upper_ug_bee,r.upper_qualifier,t) for _,r in z.iterrows()]
    z['ld50_screen']=(z.endpoint.eq('LD50') & z.unit_mapping.ne('ineligible_unit') & z.value_kind.ne('missing_or_ineligible') & z.basis_for_effect.eq('mortality'))
    # Agency identity evidence separates metiram (9006-42-2) from zineb
    # (12122-67-7). The export's polymer representation reuses the zineb key.
    # A shared key is not sufficient to equate these test substances.
    z['source_identity_status']='structure_screen_only'
    metiram_conflict=(z['CAS'].fillna('').astype(str).eq('9006-42-2') &
                      z['inchikey'].eq('AMHNZOICSMBGDH-UHFFFAOYSA-L'))
    z.loc[metiram_conflict,'source_identity_status']='ambiguous_metiram_zineb_polymer_mapping'
    z['single_structure_screen']=(z.ld50_screen & z.composition.eq('mono-constituent substance') &
        z.inchikey.fillna('').ne('') & ~z.inchi_smiles_disagree & ~metiram_conflict)
    save(z[metiram_conflict | z.inchi_smiles_disagree], 'source_identity_quarantine.csv')
    z['48h_no_material_flag_screen']=(z.single_structure_screen & z.duration_hours.eq(48) & ~z.material_review_needed & ~z.dose_basis.eq('test mat.'))
    save(z,'oft_acute_endpoints_aligned.csv')
    save(z[z.inchikey.fillna('').ne('') & ~z.overlap_inchikey],'oft_fullkey_unmatched_review.csv')
    save(z[z.conservative_unseen_candidate & z.ld50_screen],'oft_conservative_candidates.csv')
    summary=[]; thresholds=[]
    for name, subset in [('all_acute',z),('ld50_screen',z[z.ld50_screen]),('single_structure_screen',z[z.single_structure_screen]),('48h_no_material_flag_screen',z[z['48h_no_material_flag_screen']])]:
        summary.append(counts(subset,'OpenFoodTox',name,'all'))
        for route,r in subset.groupby('route'):
            summary.append(counts(r,'OpenFoodTox',name,route))
        thresholds+=threshold_tables(subset,'OpenFoodTox',name)
    info={'bee_export_rows':len(bee),'acute_export_rows':len(selected),'unique_acute_results':len(z),
        'duplicate_result_export_rows':len(selected)-len(z),'bee_endpoint_categories':bee['AdministrativeData.Endpoint'].value_counts().to_dict(),
        'acute_unit_counts':z.unit.value_counts(dropna=False).to_dict(),'acute_duration_hours':z.duration_hours.fillna('missing').value_counts().to_dict(),
        'acute_reference_structure_conflict_substances':z.loc[z.inchi_smiles_disagree,'substance_id'].nunique(),
        'novel_candidate_names':sorted(z.loc[z.conservative_unseen_candidate & z.ld50_screen,'name'].unique())}
    return z,summary,thresholds,info

def parse_value(value):
    raw=s(value).replace('≥','>=').replace('≤','<=')
    if not raw or raw in ('-','NR'):
        return None,'','missing'
    hit=re.fullmatch(r'\s*(>=|<=|>|<|=)?\s*([0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)\s*',raw)
    if not hit:
        return None,'','unparsed'
    q=hit[1] or ''
    return float(hit[2]),q,'exact' if q in ('','=') else 'censored'

def plos():
    p=pd.read_csv(OUT/'plos_identity_alignment.csv').fillna('')
    p=p.rename(columns={'Pesticide':'name'})
    # S1 supplies PCCODE but no CAS/structure. Mapping sources remain in the table.
    p['CAS']=p.get('mapping_CAS','')
    p['substance_id']=p['PCCODE']
    p=qualify(p)
    s4=pd.read_pickle(CACHE/'plos_s4.pkl')
    pc=lambda v:re.sub(r'\.0$','',str(v).strip()).zfill(6)
    refs=s4[s4.PCCODE.notna()].copy()
    refs['pc']=refs.PCCODE.map(pc)
    refs=refs.set_index('pc')
    cite=s4[s4.MRID.notna()][['MRID','Citation']].copy()
    cite['MRID']=cite.MRID.map(lambda v:re.sub(r'\.0$','',s(v)).zfill(8))
    save(cite,'plos_mrid_citations.csv')
    records=[]
    for _,row in p.iterrows():
        for route,col,mridcol in [('Contact','LD50 (μg/bee)','AAC 850.302'),('Oral','LD50 (μg/bee).1','AAO OECD 213')]:
            v,q,kind=parse_value(row[col])
            if kind=='missing':
                continue
            rec=row.to_dict()
            rec.update(route=route,endpoint='LD50',lower_ug_bee=v,lower_qualifier=q,value_kind=kind,raw_value=row[col],source_cell=('F' if route=='Contact' else 'G')+str(row.excel_row))
            hits=refs.loc[[pc(row.PCCODE)]] if pc(row.PCCODE) in refs.index else refs.iloc[:0]
            rec['MRID_mapping_status']='unique_PCCODE' if len(hits)==1 else 'unresolved'
            if len(hits)>1:
                hits=hits[hits.Pesticide.map(norm).eq(norm(row['name']))]
                rec['MRID_mapping_status']='PCCODE_and_name' if len(hits)==1 else 'ambiguous'
            rec['MRID']=s(hits.iloc[0][mridcol]) if len(hits)==1 else ''
            for t in THRESHOLDS:
                rec[f'y_{t}']=label(v,q,None,'',t)
            records.append(rec)
    z=pd.DataFrame(records)
    save(z,'plos_adult_acute_endpoints_aligned.csv')
    save(z[z.conservative_unseen_candidate],'plos_conservative_candidates.csv')
    summary=[]; threshold=[]
    for name,subset in [('all_adult_acute',z),('resolved_structure',z[z.inchikey.ne('')]),('exact_values',z[z.value_kind.eq('exact')])]:
        summary.append(counts(subset,'PLOS2022',name,'all'))
        for route,r in subset.groupby('route'):
            summary.append(counts(r,'PLOS2022',name,route))
        threshold+=threshold_tables(subset,'PLOS2022',name)
    info={'S1_chemical_rows':len(p),'S1_adult_acute_chemicals':z.substance_id.nunique(),
         'adult_acute_unresolved_names':sorted(z.loc[z.inchikey.eq(''),'name'].unique()),
         'unparsed_values':z.loc[z.value_kind.eq('unparsed'),['name','route','raw_value']].to_dict('records'),
         'S4_unique_MRID_citations':cite.MRID.nunique(),
         'adult_endpoints_missing_MRID':z.loc[z.MRID.isin(['','-']),['name','route','source_cell']].to_dict('records'),
         'novel_candidate_names':sorted(z.loc[z.conservative_unseen_candidate,'name'].unique())}
    return z,summary,threshold,info

def excluded():
    raw=pd.read_csv(OUT/'excluded_all_parsed.csv')
    u=raw[raw['Exclusion reason'].eq('Unspecifed toxicity level')].copy()
    # Upstream writes this block AFTER numeric unit conversion, but leaves the
    # original unit strings in place. Do NOT convert these numbers a second time.
    u['value_preconverted_ug_bee']=pd.to_numeric(u.observed_response_mean)
    g=u.groupby(['CAS','toxicity_type']).value_preconverted_ug_bee.agg(['min','max','median','size']).reset_index()
    for t in THRESHOLDS:
        g[f'crosses_{t}']=(g['min']<=t)&(g['max']>t)
    g['upstream_unspecified_rule']=(g['min']<11)&(g['max']>11)
    save(g,'excluded_uncertain_route_groups.csv')
    ids=pd.read_csv(OUT/'excluded_uncertain_identity_alignment.csv').fillna('')
    ac=API.set_index('CAS')
    ids['retained_same_CAS']=ids.CAS.isin(ac.index)
    ids['retained_name']=ids.CAS.map(ac.name)
    ids['retained_source']=ids.CAS.map(ac.source)
    ids['retained_route']=ids.CAS.map(ac.toxicity_type)
    ids['retained_tier']=ids.CAS.map(ac.tier)
    for t in THRESHOLDS:
        ids[f'any_route_crosses_{t}']=ids.CAS.map(g.groupby('CAS')[f'crosses_{t}'].any())
    fpgen=GetMorganGenerator(radius=2,fpSize=1024)
    afps=[fpgen.GetFingerprint(Chem.MolFromSmiles(sm)) for sm in API.canonical_smiles]
    neighbours=[]
    for _,row in ids.iterrows():
        if not row.canonical_smiles:
            continue
        fp=fpgen.GetFingerprint(Chem.MolFromSmiles(row.canonical_smiles))
        sim=np.array(DataStructs.BulkTanimotoSimilarity(fp,afps))
        # Exclude same parent and same CAS as well as exact key to prevent self-neighbours.
        same=(API.inchikey.eq(row.inchikey)|API.connectivity.eq(row.connectivity)|API.parent_connectivity.eq(row.parent_connectivity)|API.CAS.eq(row.CAS)).to_numpy()
        sim[same]=-1
        order=np.argsort(-sim,kind='stable')[:10]
        for rank,i in enumerate(order,1):
            neighbours.append({'uncertain_CAS':row.CAS,'uncertain_name':row['name'],'rank':rank,
                'neighbour_name':API.iloc[i]['name'],'neighbour_CAS':API.iloc[i].CAS,
                'neighbour_tier':API.iloc[i].tier,'tanimoto':float(sim[i]),'self_or_same_parent_removed':int(same.sum())})
    nn=pd.DataFrame(neighbours)
    save(nn,'excluded_structure_neighbours.csv')
    save(ids,'excluded_uncertain_compound_summary.csv')
    save(raw['Exclusion reason'].value_counts().rename_axis('reason').reset_index(name='records'),'excluded_reason_counts.csv')
    info={'all_excluded_records':len(raw),'blocks':raw.source_block.nunique(),'uncertain_records':len(u),
          'uncertain_CAS':u.CAS.nunique(),'uncertain_route_groups':len(g),
          'route_groups_actually_unspecified':int(g.upstream_unspecified_rule.sum()),
          'uncertain_log_route_groups_not_unspecified':int((~g.upstream_unspecified_rule).sum()),
          'retained_by_CAS':int(ids.retained_same_CAS.sum()),
          'retained_sources':ids.retained_source.value_counts().to_dict(),
          'retained_tiers':ids.retained_tier.value_counts().to_dict(),
          'unretained_CAS':ids.loc[~ids.retained_same_CAS,['CAS','name']].to_dict('records'),
          'identity_resolved':int(ids.inchikey.ne('').sum()),
          'fullkey_overlap':int(ids.overlap_inchikey.sum()),
          'nearest_nonself_tiers':nn.loc[nn['rank'].eq(1),'neighbour_tier'].value_counts().to_dict(),
          'nearest_similarity_median':float(nn.loc[nn['rank'].eq(1),'tanimoto'].median()),
          'preconverted_values_with_non_ug_unit_strings':int((~u.observed_response_unit.isin(['AI ug/org','ug/bee','ug/org'])).sum())}
    return info

def main():
    oz,os,ot,oi=oft()
    pz,ps,pt,pi=plos()
    ei=excluded()
    save(pd.DataFrame(os+ps),'overlap_summary.csv')
    save(pd.DataFrame(ot+pt),'threshold_label_feasibility.csv')
    # Shared compounds are useful for agreement checks; all three datasets may
    # still share underlying studies, so this is not independent validation.
    concord=[]
    for source,z in [('OpenFoodTox',oz[oz.single_structure_screen]),('OpenFoodTox_48h_screen',oz[oz['48h_no_material_flag_screen']]),('PLOS2022',pz)]:
        for (key,route),g in z[z.inchikey.fillna('').ne('') & z.overlap_inchikey].groupby(['inchikey','route']):
            ar=API[API.inchikey.eq(key)].iloc[0]
            if ar.toxicity_type!=route:
                continue
            rec={'source':source,'name':ar['name'],'inchikey':key,'route':route,'n_endpoint_records':len(g)}
            rec['apistox_source']=ar.source
            rec['external_names']='|'.join(sorted(set(g['name'])))
            rec['external_values_ug_bee']='|'.join(sorted(set(s(q)+s(v) for q,v in zip(g.lower_qualifier,g.lower_ug_bee))))
            rec['contains_right_censored_100']=bool((g.lower_qualifier.eq('>')&g.lower_ug_bee.eq(100)).any())
            rec['all_exact']=bool(g.value_kind.eq('exact').all())
            for t,tier in [(1,3),(11,2),(100,1)]:
                vals=set(g[f'y_{t}'])
                ext=next(iter(vals)) if len(vals)==1 and -1 not in vals else -1
                rec[f'external_y_{t}']=ext
                rec[f'apistox_y_{t}']=int(ar.tier>=tier)
                rec[f'agreement_{t}']=int(ext==rec[f'apistox_y_{t}']) if ext!=-1 else -1
            concord.append(rec)
    concord=pd.DataFrame(concord)
    save(concord,'same_route_label_concordance.csv')
    save(concord[(concord[[f'agreement_{t}' for t in THRESHOLDS]]==0).any(axis=1)],'label_disagreement_review.csv')
    cs=[]
    for source,g in concord.groupby('source'):
        fixed=g[(g[[f'agreement_{t}' for t in THRESHOLDS]]!=-1).all(axis=1)]
        for cohort,subset in [('threshold_specific',g),('all_three_determinate',fixed)]:
            for t in THRESHOLDS:
                eligible=subset[subset[f'agreement_{t}'].ne(-1)]
                cs.append({'source':source,'cohort':cohort,'threshold_ug_bee':t,'n':len(eligible),
                    'agreements':int(eligible[f'agreement_{t}'].eq(1).sum()),
                    'disagreements':int(eligible[f'agreement_{t}'].eq(0).sum()),
                    'disagreement_rate':float(eligible[f'agreement_{t}'].eq(0).mean()),
                    'discordance_with_any_right_censored_100':int((eligible[f'agreement_{t}'].eq(0)&eligible.contains_right_censored_100).sum())})
    save(pd.DataFrame(cs),'label_concordance_summary.csv')
    inventory=json.loads((Path(__file__).resolve().parents[2]/'data/external_raw/INPUT_MANIFEST.json').read_text(encoding='utf-8'))
    unchanged=all(hashlib.sha256(((Path(__file__).resolve().parents[2] / 'data/external_raw')/x['file']).read_bytes()).hexdigest()==x['sha256'] for x in inventory)
    assert unchanged
    report={'OpenFoodTox':oi,'PLOS2022':pi,'excluded':ei,
        'invariants':{'source_files_unchanged':unchanged,'no_model_fitted':True,
        'no_censored_value_imputed':True,'no_unresolved_identity_counted_as_new':True}}
    (OUT/'audit_summary.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
    print('OFT',oi)
    print('PLOS',pi)
    print('Excluded',ei)
    print(pd.DataFrame(os+ps).query("route=='all'").to_string(index=False))

if __name__=='__main__':
    main()
