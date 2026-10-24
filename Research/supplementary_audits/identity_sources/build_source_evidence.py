"""Build read-only provenance ledgers and localized Benomyl impact evidence."""
from pathlib import Path
import argparse,csv,json,hashlib,collections,statistics,datetime
from urllib.request import Request, urlopen

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
TOP=RESEARCH.parent
ROTHAMSTED_URL='https://www.era.rothamsted.ac.uk/eradoc/OCR/ResReport1977p2-57-74'
ROTHAMSTED_SHA256='270ae180f6a58f7c03133ff6b226d5b9ca2c88dd7174ece6b4d8c80e6c039ae1'
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--download-missing-source',action='store_true',
    help='Download the original Rothamsted PDF only if its packaged copy is missing; require the frozen SHA256.')
parser.add_argument('--check-source-only',action='store_true',
    help='Check the packaged Rothamsted PDF and its hash without regenerating evidence ledgers.')
args=parser.parse_args()

def resolve_rothamsted_source(download_missing=False):
    """Prefer the packaged snapshot; an optional download must match its fixed hash."""
    packaged_source=HERE/'source_documents/35461_rothamsted_official.pdf'
    if packaged_source.exists():
        content=packaged_source.read_bytes()
    elif download_missing:
        request=Request(ROTHAMSTED_URL,headers={'User-Agent':'ApisTox provenance verification'})
        with urlopen(request,timeout=60) as response:content=response.read()
    else:
        raise FileNotFoundError('Missing packaged source_documents/35461_rothamsted_official.pdf; '
            'restore it from the archive or use --download-missing-source for a hash-checked download.')
    sha=hashlib.sha256(content).hexdigest()
    if sha!=ROTHAMSTED_SHA256:
        raise ValueError('Rothamsted PDF SHA256 differs from the archived source; refusing to use or replace it.')
    if not packaged_source.exists():
        packaged_source.parent.mkdir(exist_ok=True)
        packaged_source.write_bytes(content)
    return packaged_source,sha

if args.check_source_only:
    checked_source,checked_sha=resolve_rothamsted_source(args.download_missing_source)
    print(json.dumps({'status':'passed','source':str(checked_source.relative_to(HERE)),
        'sha256':checked_sha,'ledger_rebuild_performed':False},indent=2))
    raise SystemExit(0)
def read_csv(p): return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def write_csv(p,rows):
    if not rows:return
    with p.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader();w.writerows(rows)
def write_json(p,d): p.write_text(json.dumps(d,ensure_ascii=False,indent=2),encoding='utf-8')
def safe_json(p):
    if not p.exists():return None
    try:return json.loads(p.read_text(encoding='utf-8'))
    except ValueError:return None
rawpath=RESEARCH/'data/raw/apistox.csv'
bench=read_csv(rawpath)
for i,r in enumerate(bench):r['Index']=i
groups=collections.defaultdict(list)
for r in bench:groups[r['CID']].append(r)
known=[r for g in groups.values() if len(g)>1 for r in g]
identity={int(r['Index']):r for r in read_csv(HERE/'cid_structure_ledger_1035.csv')}
member={}
for split in ['random','maxmin','time']:
    for part in ['train','test']:
        for r in read_csv(RESEARCH/f'data/official_splits/{split}_{part}.csv'):
            member[(split,r['CAS'],r['SMILES'])]=part

def find_sections(items,heading):
    out=[]
    for s in items:
        if s.get('TOCHeading')==heading:out.append(s)
        out.extend(find_sections(s.get('Section',[]),heading))
    return out

def create_date(cid):
    j=safe_json(HERE/'responses'/f'pugview_{cid}.body') or {}
    sections=find_sections(j.get('Record',{}).get('Section',[]),'Create Date')
    try:return sections[0]['Information'][0]['Value']['DateISO8601'][0]
    except (KeyError,IndexError):return ''

known_rows=[]
for row in known:
    i=row['Index'];base=identity[i];cas=row['CAS'].strip()
    candidates=(safe_json(HERE/'responses'/f'cas_{cas}.body') or {}).get('PropertyTable',{}).get('Properties',[])
    all_candidates=';'.join(str(r['CID']) for r in candidates)
    full_candidates=[r for r in candidates if r['InChIKey']==base['benchmark_inchikey']]
    conn_candidates=[r for r in candidates if r['InChIKey'].split('-')[0]==base['benchmark_inchikey'].split('-')[0]]
    if full_candidates: candidate_status='CAS_lookup_and_frozen_SMILES_full_InChIKey_agree'
    elif conn_candidates: candidate_status='CAS_lookup_connectivity_agrees_full_key_differs'
    elif cas=='9006-42-2':candidate_status='non_discrete_reference_substance_no_current_CID'
    else:candidate_status='no_agreeing_current_CAS_candidate'
    if i in [413,672,832]:assessment='current_CID_refers_to_different_named_structure;CAS_candidate_matches_frozen_SMILES'
    elif i==290:assessment='cation_vs_chloride_salt;CAS_candidate_matches_frozen_cation'
    elif i==323:assessment='metiram_is_distinct_from_zineb;shared_fragment_does_not_identify_material;no_discrete_CID_assignment'
    elif i==824:assessment='connectivity_agrees;benchmark_unspecified_oxime_geometry_vs_current_E_isomer'
    else:assessment='current_CID_full_InChIKey_matches_frozen_SMILES'
    cid=int(row['CID']); lit_meta=safe_json(HERE/'responses'/f'literature_{cid}.json') or {}
    full_json=safe_json(HERE/'responses'/f'pugview_{cid}.body') or {}
    ag_sections=find_sections(full_json.get('Record',{}).get('Section',[]),'Agrochemical Information')
    rec=dict(row,comparison=base['comparison'],benchmark_inchikey=base['benchmark_inchikey'],
        stored_cid_current_inchikey=base['current_cid_inchikey'],current_CAS_candidate_CIDs=all_candidates,
        full_key_agreeing_candidate_CIDs=';'.join(str(r['CID']) for r in full_candidates),
        connectivity_agreeing_candidate_CIDs=';'.join(str(r['CID']) for r in conn_candidates),
        candidate_evidence_status=candidate_status,assessment=assessment,
        current_CAS_response='responses/cas_'+cas+'.body',
        current_CAS_response_sha256=hashlib.sha256((HERE/'responses'/('cas_'+cas+'.body')).read_bytes()).hexdigest(),
        stored_CID_current_create_date=create_date(cid),
        current_candidate_create_dates=';'.join(str(r['CID'])+':'+create_date(r['CID']) for r in candidates),
        literature_endpoint_http_status=lit_meta.get('status','not_queried'),
        current_earliest_literature_date='unresolved;old_SDQ_endpoint_did_not_return_literature_JSON',
        frozen_year_origin='PubChem_literature_or_create_date;per_row_historical_response_not_preserved',
        frozen_year_error_established=False,year_replacement_authorized_by_evidence=False,
        time_membership_error_established=False,
        use_flag_source='CID_PubChem_Agrochemical_Information' if row['source']=='ECOTOX' else 'PPDB_BPDB_HTML_text',
        current_Agrochemical_Information_heading_present=bool(ag_sections),
        frozen_use_error_from_CID_established=False,
        use_effect_note='CID_dependency_exists_but_current_CID_matches_this_ECOTOX_row' if row['source']=='ECOTOX' else 'Uses_came_from_provider_HTML;CID_mismatch_does_not_directly_change_use_flags',
        manual_expert_dual_adjudication=False,frozen_record_changed=False)
    for split in ['random','maxmin','time']:rec[split+'_membership']=member[(split,row['CAS'],row['SMILES'])]
    known_rows.append(rec)
write_csv(HERE/'known_5_CID_groups_10_rows_impact_ledger.csv',known_rows)

candidate_rows=[]
for flagged in read_csv(HERE/'current_identifier_flagged.csv'):
    cas=flagged['CAS'].strip()
    content=safe_json(HERE/'responses'/f'cas_{cas}.body') or {}
    candidates=content.get('PropertyTable',{}).get('Properties',[])
    full=[c for c in candidates if c['InChIKey']==flagged['benchmark_inchikey']]
    connectivity=[c for c in candidates if c['InChIKey'].split('-')[0]==flagged['benchmark_inchikey'].split('-')[0]]
    rec=dict(flagged,current_CAS_candidate_CIDs=';'.join(str(c['CID']) for c in candidates),
        full_key_agreeing_candidate_CIDs=';'.join(str(c['CID']) for c in full),
        connectivity_agreeing_candidate_CIDs=';'.join(str(c['CID']) for c in connectivity),
        current_CAS_response='responses/cas_'+cas+'.body',
        current_CAS_response_sha256=hashlib.sha256((HERE/'responses'/('cas_'+cas+'.body')).read_bytes()).hexdigest(),
        candidate_status=('current_CAS_candidate_matches_frozen_SMILES' if full else
            'current_CAS_candidate_connectivity_matches_frozen_SMILES' if connectivity else
            'CAS_returned_structures_require_further_reconciliation' if candidates else
            'CAS_lookup_returned_no_discrete_CID'),
        approved_historical_replacement=False,expert_dual_adjudication=False)
    candidate_rows.append(rec)
write_csv(HERE/'flagged_33_CAS_candidate_ledger.csv',candidate_rows)

# Establish which logged observations belonged to actually rejected route groups.
factor={'AI ug/org':1,'AI ug/org/d':1,'ug/bee':1,'ug/org':1,'ug/org/d':1,
        'AI ng/org':.001,'AI ng/org/d':.001,'ng/org':.001,'AI mg/org':1000,
        'mg/bee':1000,'mg/org':1000,'pg/org':.000001}
route={'Diet, unspecified':'Oral','Drinking water':'Oral','Food':'Oral','Dermal':'Contact',
       'Direct application':'Contact','Topical, general':'Contact',
       'Multiple routes between application groups':'Other','Oral via capsule':'Other',
       'Spray, unspecified':'Other','Environmental, unspecified':'Other'}
ecopath=RESEARCH/'output/external_data_audit/cache/upstream_ecotox.csv'
normalized=[];egroups=collections.defaultdict(list)
with ecopath.open(encoding='utf-8-sig') as f:
    for idx,original in enumerate(csv.DictReader(f,delimiter='|')):
        r={k.strip():v.strip() for k,v in original.items()}
        u=r['Observed Response Units'];v=r['Observed Response Mean'];rt=route.get(r['Exposure Type'])
        if u not in factor or v=='NR' or rt is None:continue
        c=r['CAS Number'];cas=f'{c[:-3]}-{c[-3:-1]}-{c[-1]}'
        rec={'raw_row_0based':idx,'CAS':cas,'route':rt,'value_raw':v,'unit_raw':u,
             'value_converted':float(v.replace('/',''))*factor[u],
             'operator':r['Observed Response Mean Op']}
        normalized.append(rec);egroups[cas,rt].append(rec)
uncertain={key for key,value in egroups.items()
           if min(r['value_converted'] for r in value)<11<max(r['value_converted'] for r in value)}
ucas={key[0] for key in uncertain}
logs=[]
for r in normalized:
    if r['CAS'] not in ucas:continue
    rec=dict(r,actually_rejected_route_group=(r['CAS'],r['route']) in uncertain,
             included_in_upstream_exclusion_log=True)
    rec['log_status']='rejected_route_group' if rec['actually_rejected_route_group'] else 'retained_route_group_logged_due_to_shared_CAS'
    logs.append(rec)
write_csv(HERE/'upstream_exclusion_log_745_row_audit.csv',logs)
summary={'log_rows':len(logs),'log_CAS':len(ucas),'actual_rejected_groups':len(uncertain),
         'actual_rejected_records':sum(r['actually_rejected_route_group'] for r in logs),
         'retained_group_rows_also_logged':sum(not r['actually_rejected_route_group'] for r in logs),
         'retained_groups_also_logged':len({(r['CAS'],r['route']) for r in logs if not r['actually_rejected_route_group']}),
         'ecotox_snapshot_sha256':hashlib.sha256(ecopath.read_bytes()).hexdigest(),
         'evidence':'upstream/dataset_creation_ecotox.py:175-186;log masks by CAS while filtering masks route statistics',
         'model_impact':'The log defect does not change the route-level filtering executed by the original code.'}
assert summary['log_rows']==745 and summary['actual_rejected_records']==536 and summary['retained_group_rows_also_logged']==209
write_json(HERE/'exclusion_log_summary.json',summary)

# Benomyl source discrepancy is a proposed interpretation, not an adjudicated correction.
ben=[dict(r) for r in normalized if r['CAS']=='17804-35-2']
for r in ben:
    r['upstream_eligible_route']=(r['CAS'],r['route']) not in uncertain
    r['source_discrepancy_target']=r['raw_row_0based']==1251
    r['proposed_operator_if_original_table_mapping_accepted']='>' if r['source_discrepancy_target'] else r['operator']
    r['acceptance_status']='not_expert_adjudicated;original_trial_identity_not_fully_established'
write_csv(HERE/'benomyl_record_level_impact_ledger.csv',ben)
assert len(ben)==7 and sum(r['upstream_eligible_route'] for r in ben)==5
def interval(r,proposed):
    v=r['value_converted'];op=r['proposed_operator_if_original_table_mapping_accepted'] if proposed else r['operator']
    if op=='>':return v,float('inf'),True
    if op=='>=':return v,float('inf'),False
    if op in ('<','<='):return 0.,v,False
    return v,v,False
def classify_interval(lo,hi,open_lo,t):
    if hi<=t:return 1
    if lo>t or (lo==t and open_lo):return 0
    return -1
def benomyl_label(t,rule,proposed):
    by_route=collections.defaultdict(list)
    for r in ben:
        if r['upstream_eligible_route']:by_route[r['route']].append(r)
    route_labels=[]
    for rr in by_route.values():
        intervals=[interval(r,proposed) for r in rr]
        if rule=='C':
            labels={classify_interval(lo,hi,op,t) for lo,hi,op in intervals}
            route_labels.append(next(iter(labels)) if len(labels)==1 and -1 not in labels else -1)
        else:
            lo=statistics.median(x[0] for x in intervals);hi=statistics.median(x[1] for x in intervals)
            # No Benomyl median-bound equality with the three thresholds occurs.
            # This local check does not claim to replace the study's general interval-median implementation.
            assert lo!=t and hi!=t
            route_labels.append(classify_interval(lo,hi,False,t))
    return 1 if 1 in route_labels else 0 if all(x==0 for x in route_labels) else -1
effects=[]
for threshold,rule,original in [(100,'D',1),(11,'C',1),(1,'D',0)]:
    old=benomyl_label(threshold,rule,False);new=benomyl_label(threshold,rule,True)
    effects.append({'threshold':threshold,'primary_rule':rule,'historical_eligible_route_label':old,
        'proposed_eligible_route_label':new,'keep_old':old if old!=-1 else original,
        'keep_proposed':new if new!=-1 else original})
assert [(r['historical_eligible_route_label'],r['proposed_eligible_route_label']) for r in effects]==[(1,-1),(1,-1),(-1,0)]
for r in effects:
    r['CAS']='17804-35-2';r['Index']=111;r['exclude_membership_old']='excluded_due_to_unresolved_at_1'
    r['exclude_membership_proposed']='excluded_due_to_unresolved_at_11_and_100'
    r['status']='conditional_rule_propagation_only;not_approved_label_replacement'
write_csv(HERE/'benomyl_primary_rule_conditional_effects.csv',effects)
packaged_source,sha=resolve_rothamsted_source(args.download_missing_source)
write_json(HERE/'source_documents/35461_provenance.json',{
    'official_url':'https://www.era.rothamsted.ac.uk/eradoc/OCR/ResReport1977p2-57-74',
    'doi':'10.23637/ERADOC-1-16','sha256':sha,'table':'Table 2','PDF_page':3,'printed_page':56,
    'visual_check':'Benomyl contact >10 and oral >10 ug/bee; comparison with historical ECOTOX oral <10',
    'reused_local_original':str(packaged_source.relative_to(HERE)),'matches_latest_trace_expected_sha256':True,
    'adjudication':'Document-level discrepancy; one-to-one original trial identity not fully established; no named two-expert independent adjudication claimed.'})
write_json(HERE/'source_evidence_summary.json',{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'known_CID_groups':5,'known_rows':10,'known_noncorresponding_distinct_current_structures':[413,672,832],
    'salt_form_case':[290],'non_discrete_metiram_case':[323],'oxime_geometry_case':[824],
    'year_errors_proven':0,'time_membership_changes_proven':0,'use_errors_proven_from_CID':0,
    'zero_proven_errors_note':'Zero established downstream errors means unresolved evidence, not demonstrated correctness.',
    'frozen_inputs_changed':False,'model_fitting_performed_by_this_script':False,
    'exclusion_log':summary})
print(json.dumps({'known_rows':len(known_rows),'exclusion_log':summary},indent=2))
