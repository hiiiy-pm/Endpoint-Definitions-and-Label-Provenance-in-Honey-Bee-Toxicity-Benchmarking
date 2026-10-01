"""Round 2: rank-based metrics, censoring, identity exclusion and figure checks.
No production evaluator or sklearn metric functions are imported.
"""
from pathlib import Path
import re,json
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
QA=ROOT/'verification/results'; QA.mkdir(parents=True,exist_ok=True)
checks={};errors=[]
def check(k,v):
 checks[k]=bool(v)
 if not v: errors.append(k)
def equal(a,b):return bool(np.isclose(a,b,atol=1e-12,rtol=0,equal_nan=True))
def ranking(y,z):
 y=np.asarray(y);z=np.asarray(z);pos=z[y==1];neg=z[y==0]
 if len(pos)==0 or len(neg)==0:return np.nan,np.nan
 auc=((pos[:,None]>neg).sum()+.5*(pos[:,None]==neg).sum())/(len(pos)*len(neg))
 ap=np.mean([np.sum((z>=p)&(y==1))/np.sum(z>=p) for p in pos])
 return float(auc),float(ap)
# Every primary metric is independently checked by positive-negative rank pairs.
p=pd.read_csv(ROOT/'results/primary/02_test_predictions.csv')
for a in pd.read_csv(ROOT/'results/primary/01_representation_performance.csv').itertuples():
 d=p[(p.Split==a.Split)&(p.Representation==a.Representation)&(p.Threshold==a.Threshold)]
 auc,ap=ranking(d.Y,d.Score)
 check('rank_primary:'+a.Split+'/'+a.Representation+'/'+str(a.Threshold),equal(auc,a.AUROC) and equal(ap,a.AUPRC))
# Cross-source metrics and every stored bootstrap interval and count.
f=ROOT/'output/external_evaluation'
scores=pd.read_csv(f/'frozen_scores_and_labels.csv');scores=scores[scores.Representation=='ECFP']
meta=pd.read_csv(f/'frozen_compound_members.csv');metrics=pd.read_csv(f/'frozen_metrics_ECFP.csv');contrasts=pd.read_csv(f/'frozen_contrasts_ECFP.csv')
npz=np.load(f/'frozen_ECFP_paired_bootstrap_distributions.npz');computed={}
for a in metrics.itertuples():
 key=(a.source,a.split,a.cohort,a.labels,a.threshold)
 q=meta[(meta.source==a.source)&(meta['split']==a.split)]
 if 'original_contact' in a.cohort:q=q[q.original_route=='Contact']
 if 'quality_screen' in a.cohort:q=q[q.passes_plos_quality_screen.eq(True)]
 d=scores[(scores.source==a.source)&(scores['split']==a.split)&(scores.threshold==a.threshold)&scores.inchikey.isin(q.inchikey)]
 auc,ap=ranking(d[a.labels+'_label'],d.Score)
 computed[key]=auc
 tag='/'.join(map(str,key))+'/'+a.metric
 check('rank_source:'+tag,equal(auc if a.metric=='AUROC' else ap,a.estimate))
 npos=int(d[a.labels+'_label'].sum());check('sparse_flags:'+tag,bool(a.exploratory_sparse)==(min(npos,len(d)-npos)<5))
 payload=npz['__'.join(map(str,[a.source,a.split,a.cohort,a.labels,a.threshold,a.metric]))];v=payload[np.isfinite(payload)]
 lo,hi=np.percentile(v,[2.5,97.5]) if len(v) else (np.nan,np.nan)
 check('source_bootstrap:'+tag,len(v)==a.bootstrap_valid and len(payload)-len(v)==a.bootstrap_degenerate and equal(lo,a.ci_low) and equal(hi,a.ci_high))
for i,a in enumerate(contrasts.itertuples()):
 t=11 if '11-' in a.comparison else 1
 def delta(l):return computed[(a.source,a.split,a.cohort,l,t)]-computed[(a.source,a.split,a.cohort,l,100)]
 expected=delta(a.labels) if a.labels in ('original','external') else delta('external')-delta('original')
 check('source_contrast:'+str(i),equal(expected,a.estimate))
# Unseen-identity case study: interval semantics and exclusion from the benchmark.
f=ROOT/'output/unseen_identity_case_study'
e=pd.read_csv(f/'data_curation/frozen_scoreable_48h_endpoints.csv')
ref=pd.read_csv(f/'model_restore/original_structure_reference.csv')
def classify(q,v,t):
 q='' if pd.isna(q) else str(q).strip()
 if q in ('','=', 'exact'):return int(v<=t)
 if q=='>':return 0 if v>=t else np.nan
 if q=='>=':return 0 if v>t else np.nan
 if q in ('<','<='):return 1 if v<=t else np.nan
 raise ValueError(q)
for a in e.itertuples():
 for t in [1,11,100]:check('censor_label:'+a.endpoint_id+'/'+str(t),equal(classify(a.qualifier,a.LD50_value,t),getattr(a,'label_le_'+str(t))))
 for key in ['inchikey','connectivity','parent_connectivity']:check('unseen_identity:'+a.endpoint_id+'/'+key,getattr(a,key) not in set(ref[key].dropna()))
 check('endpoint_48h:'+a.endpoint_id,a.observation_h==48 and a.species=='Apis mellifera')
check('six_unique_identities_eleven_endpoints',e.compound_id.nunique()==6 and len(e)==11 and e.route.value_counts().to_dict()=={'oral':6,'contact':5})
# Recompute all 72 case-study metrics and confusion matrices on the stated populations.
case=pd.read_csv(f/'evaluation/new_compound_case_predictions.csv');perf=pd.read_csv(f/'evaluation/pilot_performance.csv')
for i,a in enumerate(perf.itertuples()):
 q=e[e.route==a.route]
 if a.quality_screen=='explicit_adult_only':q=q[q.life_stage_evidence=='explicit_in_source']
 if a.cohort=='common_three_thresholds':q=q[q[['label_le_1','label_le_11','label_le_100']].notna().all(axis=1)]
 q=q[q['label_le_'+str(a.threshold)].notna()]
 d=case[(case.endpoint_id.isin(q.endpoint_id))&(case['split']==a.split)&(case.threshold==a.threshold)]
 y=d.experimental_label.to_numpy();z=d.score.to_numpy();binary=(z>0).astype(int)
 auc,ap=ranking(y,z)
 conf=[int(((binary==b)&(y==v)).sum()) for b,v in [(1,1),(0,0),(1,0),(0,1)]]
 check('pilot_metric:'+str(i),equal(auc,a.AUROC) and equal(ap,a.average_precision))
 check('pilot_population:'+str(i),len(d)==a.n_compounds and y.sum()==a.n_positive)
 check('pilot_confusion:'+str(i),conf==[a.TP,a.TN,a.FP,a.FN])
 check('pilot_all_sparse:'+str(i),a.exploratory_sparse and a.single_class==(len(set(y))<2))
# Restored held-out scores and panel counts.
restore=pd.read_csv(f/'model_restore/restored_vs_original_test_scores.csv')
num=[c for c in restore if 'score' in c.lower()];print('Restoration score columns:',num)
if {'original_score','restored_score'}<=set(restore):
 check('1863_restored_scores',len(restore)==1863 and np.max(np.abs(restore.original_score-restore.restored_score))<=3e-16)
else:
 vals=restore[num].select_dtypes(include='number')
 check('1863_restored_scores',len(restore)==1863 and vals.shape[1]>=2 and all(np.max(np.abs(vals.iloc[:,0]-vals.iloc[:,j]))<=3e-16 for j in range(1,vals.shape[1])))
allpred=pd.concat([pd.read_csv(f/f'predictions/{d}/new_compound_predictions.csv') for d in ['initial8','added2']])
check('all10_candidates_90scores',len(allpred)==90 and allpred.compound_id.nunique()==10)
check('external_inputs_are_structure_only',all(set(pd.read_csv(f/f'predictions/{d}/structure_only_input.csv').columns)<= {'compound_id','SMILES','name','inchikey','CAS','source'} for d in ['initial8','added2']))
check('default_contact_11_all_missed',case[(case.route=='contact')&(case.threshold==11)&(case.experimental_label==1)].predicted_label.eq(0).all())
check('default_oral_severe_all_missed',case[(case.route=='oral')&case.threshold.isin([1,11])&(case.experimental_label==1)].predicted_label.eq(0).all())
# Supplementary tables S1-S20 appear once each, in citation order.
ordered=(ROOT/'output/latex_tables/supplementary_tables_ordered.tex').read_text(encoding='utf-8')
found=re.findall(r'Supplementary Table S(\d+[ab]?)\.',ordered)
expected=[str(i) for i in range(1,18)]+['18a','18b','19a','19b','20']
check('supp_tables_complete_and_ordered',found==expected)
# Figure source data and figure QA reports.
csv=pd.read_csv(ROOT/'source_data/figure_external_reassessment.csv')
check('all_24_paired_figure_rows',len(csv)==48 and csv.groupby(['block','split','comparison']).size().eq(2).all())
QA_FIG=ROOT/'output/figure_qa'
for f in sorted((ROOT/'figures').glob('Fig*.pdf')):
 qa=json.loads((QA_FIG/(f.stem+'.collisions.json')).read_text(encoding='utf-8'))
 check('figure_collision_gate:'+f.stem,qa['verdict']=='PASS' and qa['summary']['fail']==0)
 al=QA_FIG/(f.stem+'.alignment.json')
 if al.exists():
  a=json.loads(al.read_text(encoding='utf-8'));check('figure_alignment_gate:'+f.stem,a['verdict'] in ('PASS','NOT APPLICABLE') and a['summary']['fail']==0)
# Independent rank-pair recomputation of the Tier1-excluded refit AUROCs.
rf=pd.read_csv(ROOT/'results/tier_boundary/02b_tier1_refit_test_predictions.csv')
ex=pd.read_csv(ROOT/'results/tier_boundary/02_tier1_exclusion_auroc.csv');ex=ex[(ex.Bootstrap=='ScaffoldCluster')&(ex.Quantity=='Aclean_noT1_train_eval')]
for a in ex.itertuples():
 r=rf[(rf.Split==a.Split)&(rf.Representation==a.Representation)];auc,_=ranking(r.Y,r.Score)
 check('rank_tier1_refit:'+a.Split+'/'+a.Representation,equal(auc,a.AUROC))
report={'round':2,'status':'PASS_WITH_LIMITATIONS' if not errors else 'FAIL','checks_total':len(checks),'checks_passed':sum(checks.values()),'errors':errors,'checks':checks,'method':'independent positive-negative pair rankings and precision-at-score-cutoff; explicit interval semantics; all case-study confusion matrices; stored bootstrap quantiles; supplementary table order and figure QA reports','new_model_fits':0,'limits':['Selected sparse unseen-identity case study','Retrospective source reassessment, possible shared studies','Qualifier tracing limited to ECOTOX-derived records','Tier1 exclusion also removes adjacent-tier comparisons']}
(QA/'round2.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2,ensure_ascii=False));assert not errors,errors
