"""Round 1: independent metric, table-row and preservation checks."""
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score
import json
ROOT=Path(__file__).resolve().parents[1]
QA=ROOT/'verification/results';QA.mkdir(parents=True,exist_ok=True)
checks={};errors=[]
def check(name,ok):
 checks[name]=bool(ok)
 if not ok:errors.append(name)
def metric(y,z):return (roc_auc_score(y,z),average_precision_score(y,z)) if len(set(y))==2 else (np.nan,np.nan)
def close(a,b):return bool(np.isclose(a,b,rtol=0,atol=1e-12,equal_nan=True))
# Primary scores, all 45 settings.
p=pd.read_csv(ROOT/'results/primary/02_test_predictions.csv')
t=pd.read_csv(ROOT/'results/primary/01_representation_performance.csv')
for a in t.itertuples():
 d=p[(p.Split==a.Split)&(p.Representation==a.Representation)&(p.Threshold==a.Threshold)]
 auc,ap=metric(d.Y,d.Score)
 check(f'primary_auc:{a.Split}/{a.Representation}/{a.Threshold}',close(auc,a.AUROC))
 check(f'primary_ap:{a.Split}/{a.Representation}/{a.Threshold}',close(ap,a.AUPRC))
 check(f'primary_n:{a.Split}/{a.Representation}/{a.Threshold}',len(d)==207 and d.Y.sum()==a.n_pos)
# Every external metric on saved scores, independently of the production evaluator.
f=ROOT/'output/external_evaluation'
members=pd.read_csv(f/'frozen_compound_members.csv')
m=pd.read_csv(f/'frozen_metrics_ECFP.csv');scores=pd.read_csv(f/'frozen_scores_and_labels.csv');scores=scores[scores.Representation=='ECFP']
for i,a in enumerate(m.itertuples()):
 q=members[(members.source==a.source)&(members['split']==a.split)]
 if a.cohort in ('original_contact_only','original_contact_quality_screen'):q=q[q.original_route.eq('Contact')]
 if a.cohort in ('quality_screen','original_contact_quality_screen'):q=q[q.passes_plos_quality_screen]
 d=scores[(scores.source==a.source)&(scores['split']==a.split)&(scores.threshold==a.threshold)&scores.inchikey.isin(q.inchikey)]
 y=d[a.labels+'_label'];auc,ap=metric(y,d.Score)
 check(f'cross_source_metric:{i}',close(auc if a.metric=='AUROC' else ap,a.estimate))
 check(f'cross_source_counts:{i}',len(d)==a.n_compounds and y.sum()==a.n_pos and len(d)-y.sum()==a.n_neg)
# Preserved scores equal the held-out scores; row identities are excluded from training.
s=scores.merge(p[p.Representation=='ECFP'].assign(split=lambda x:x.Split.str.lower()).rename(columns={'Threshold':'threshold'}),on=['split','Representation','threshold','Index'],suffixes=('_external','_primary'),validate='many_to_one')
check('saved_scores_unchanged',np.allclose(s.Score_external,s.Score_primary,rtol=0,atol=1e-14))
for sp,d in scores.groupby('split'):
 train=pd.read_csv(ROOT/f'data/official_splits/{sp}_train.csv')
 ref=pd.read_csv(ROOT/'results/cache/processed_data.csv')
 train_ids=set(ref.index[ref.SMILES.isin(train.SMILES)])
 check('training_ids_absent:'+sp,not bool(set(d.Index)&train_ids))
# Route-matched concordance reconstructed from row-level records.
f=ROOT/'output/external_data_audit';d=pd.read_csv(f/'same_route_label_concordance.csv');c=pd.read_csv(f/'label_concordance_summary.csv')
for a in c.itertuples():
 q=d[d.source.eq(a.source)]
 if a.cohort=='all_three_determinate':q=q[(q[['external_y_1','external_y_11','external_y_100']]>=0).all(axis=1)]
 else:q=q[q[f'external_y_{a.threshold_ug_bee}']>=0]
 check(f'concordance_n:{a.source}/{a.cohort}/{a.threshold_ug_bee}',len(q)==a.n)
 check(f'concordance_disagreement:{a.source}/{a.cohort}/{a.threshold_ug_bee}',(q[f'external_y_{a.threshold_ug_bee}']!=q[f'apistox_y_{a.threshold_ug_bee}']).sum()==a.disagreements)
# Main-text Tables 3-5 are regenerated from saved outputs and compared with the written rows.
import importlib.util
_spec=importlib.util.spec_from_file_location('main_tables',ROOT/'scripts/make_main_tables.py');mt=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(mt)
for tname,rows in [('table3_tier1_exclusion',mt.table3_rows()),('table4_concordance_by_tier',mt.table4_rows()),('table5_external_auroc',mt.table5_rows())]:
 saved=(ROOT/'output/main_tables'/(tname+'.tex')).read_text(encoding='utf-8').splitlines()
 check('main_table_rows:'+tname,rows==saved)
# Compact descriptor table (Supplementary Table S8) against the saved descriptor results.
si=(ROOT/'output/latex_tables/supp_tables_external.tex').read_text(encoding='utf-8')
q=pd.read_csv(ROOT/'results/supplementary/02_topk_descriptor_performance.csv')
for (sp,k),d in q.groupby(['Split','TopK']):
 d=d.set_index('Threshold');row=f'{sp} & {k} & '+' & '.join(f'{d.loc[t,"TestAUROC"]:.3f}' for t in [100,11,1])
 check('supp_descriptor_table:'+sp+'/'+str(k),row in si)
# Tier-boundary analyses: tier-pair AUROC, the Eq. (4) identity, Tier1 exclusion and concordance.
proc=pd.read_csv(ROOT/'results/cache/processed_data.csv');tier=proc.Tier.to_numpy()
tp=pd.read_csv(ROOT/'results/tier_boundary/01_tier_pair_auroc.csv')
for a in tp.itertuples():
 d=p[(p.Split==a.Split)&(p.Representation==a.Representation)&(p.Threshold==a.TrainedFor)];t=tier[d.Index.to_numpy()]
 k=(t==a.HigherTier)|(t==a.LowerTier)
 check(f'tier_pair:{a.Split}/{a.Representation}/{a.TrainedFor}/{a.Pair}',close(roc_auc_score((t[k]==a.HigherTier).astype(int),d.Score.to_numpy()[k]),a.AUROC))
for (sp,rep,trained),g in tp.groupby(['Split','Representation','TrainedFor']):
 d=p[(p.Split==sp)&(p.Representation==rep)&(p.Threshold==trained)];t=tier[d.Index.to_numpy()]
 for thr,kb in [(100,1),(11,2),(1,3)]:
  check(f'tier_pair_identity:{sp}/{rep}/{trained}/{thr}',abs((g[f'Weight_le{thr}']*g.AUROC).sum()-roc_auc_score((t>=kb).astype(int),d.Score))<1e-9)
ex=pd.read_csv(ROOT/'results/tier_boundary/02_tier1_exclusion_auroc.csv');ex=ex[ex.Bootstrap=='ScaffoldCluster']
rf=pd.read_csv(ROOT/'results/tier_boundary/02b_tier1_refit_test_predictions.csv')
dc=pd.read_csv(ROOT/'results/tier_boundary/03_tier1_exclusion_decomposition.csv');dc=dc[dc.Bootstrap=='ScaffoldCluster']
for (sp,rep),g in ex.groupby(['Split','Representation']):
 q=g.set_index('Quantity').AUROC;r=rf[(rf.Split==sp)&(rf.Representation==rep)]
 d=p[(p.Split==sp)&(p.Representation==rep)&(p.Threshold==100)];t=tier[d.Index.to_numpy()];k=t!=1
 check(f'tier1_refit_auc:{sp}/{rep}',close(roc_auc_score(r.Y,r.Score),q['Aclean_noT1_train_eval']) and set(r.Tier)=={0,2,3})
 check(f'tier1_eval_auc:{sp}/{rep}',close(roc_auc_score((t[k]>=2).astype(int),d.Score.to_numpy()[k]),q['A100_noT1_eval']))
 c=dc[(dc.Split==sp)&(dc.Representation==rep)].set_index('Contrast').Estimate
 check(f'tier1_decomposition_identity:{sp}/{rep}',abs(c['Gap_11_minus_100']-c['Tier1_total_component']-c['Residual_11_minus_clean'])<1e-12 and abs(c['Tier1_total_component']-c['Evaluation_component']-c['Training_component'])<1e-12)
bt=pd.read_csv(ROOT/'results/tier_boundary/04_concordance_by_tier.csv')
cc=pd.read_csv(ROOT/'output/external_data_audit/same_route_label_concordance.csv')
cc=cc[(cc[['external_y_1','external_y_11','external_y_100']]>=0).all(axis=1)]
src_name=cc.source.map({'OpenFoodTox':'OFT','OpenFoodTox_48h_screen':'OFT 48 h','PLOS2022':'EPA'});ctier=cc.apistox_y_100+cc.apistox_y_11+cc.apistox_y_1
for a in bt.itertuples():
 q=cc[(src_name==a.Source)&(ctier==a.Tier)]
 check(f'concordance_tier:{a.Source}/{a.Tier}',len(q)==a.n and int((q.agreement_100==1).sum())==a.Agree_100 and int(((q.agreement_100==0)&q.contains_right_censored_100).sum())==a.Disagree_with_right_censored)
qs=json.loads((ROOT/'results/tier_boundary/QUALIFIER_PROPAGATION_SUMMARY.json').read_text(encoding='utf-8'))
check('qualifier_replication_all_441',qs['benchmark_ecotox_compounds']==441 and qs['replicated_compounds']==441 and qs['replication_label_and_level_match']==441)
for r in qs['by_tier']:
 check(f"qualifier_categories_sum:{r['Tier']}",r['interval_label_100_negative']+r['interval_label_100_unresolved']+r['interval_label_100_positive']==r['n_ecotox_compounds'])
report={'round':1,'status':'PASS' if not errors else 'FAIL','checks_passed':sum(checks.values()),'checks_total':len(checks),'errors':errors,'checks':checks,'method':'sklearn point metrics from saved prediction rows; row-level concordance; regenerated table rows','new_model_fits':0,'repository_relative_inputs':True}
(QA/'round1.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ('checks','errors')},indent=2,ensure_ascii=False))
assert not errors,errors
