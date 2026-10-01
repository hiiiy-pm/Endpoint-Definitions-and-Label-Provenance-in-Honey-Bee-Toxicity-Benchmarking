"""First-round independent checks for the new external evaluation code."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
from evaluate_new_compound_pilot import weighted_auc_many, truth, exact_binomial

OUT=Path(__file__).resolve().parent
E=OUT/'evaluation'
checks=[]
def check(name,condition):
    if not bool(condition):
        raise AssertionError(name)
    checks.append(name)

rng=np.random.default_rng(914)
max_error=0.
for i in range(100):
    y=rng.integers(0,2,size=8)
    s=rng.integers(-2,3,size=8).astype(float)
    w=rng.multinomial(8,np.full(8,1/8),size=10)
    values=weighted_auc_many(y,s,w)
    for j,ww in enumerate(w):
        if ww[y==1].sum() and ww[y==0].sum():
            expected=roc_auc_score(y,s,sample_weight=ww)
            max_error=max(max_error,abs(values[j]-expected))
            check(f'weighted_auc_{i}_{j}',abs(values[j]-expected)<1e-12)
        else:
            check(f'weighted_single_class_{i}_{j}',np.isnan(values[j]))
for value,q,t,expected in [(100,'>',100,0),(100,'>=',100,np.nan),(98.8,'>',100,np.nan),
                          (2.68,'=',1,0),(2.68,'=',11,1),(.29,'=',1,1)]:
    got=truth(value,q,t)
    check(f'interval_{value}_{q}_{t}',(np.isnan(got) and np.isnan(expected)) or got==expected)
lo,hi=exact_binomial(0,1)
check('exact_ci_zero_of_one',lo==0 and abs(hi-.975)<1e-12)
lo,hi=exact_binomial(1,1)
check('exact_ci_one_of_one',abs(lo-.025)<1e-12 and hi==1)
check('empty_exact_ci',all(np.isnan(exact_binomial(0,0))))

metrics=pd.read_csv(E/'pilot_performance.csv')
cases=pd.read_csv(E/'new_compound_case_predictions.csv')
check('metric_grid',len(metrics)==72)
check('case_grid',len(cases)==99)
check('class_totals',metrics.n_compounds.eq(metrics.n_positive+metrics.n_negative).all())
check('confusion_totals',(metrics[['TP','TN','FP','FN']].sum(axis=1)==metrics.n_compounds).all())
check('bootstrap_totals',(metrics.bootstrap_valid+metrics.bootstrap_degenerate).eq(5000).all())
check('single_class_auc_na',metrics.loc[metrics.single_class,'AUROC'].isna().all())
check('single_class_ap_na',metrics.loc[metrics.single_class,'average_precision'].isna().all())
check('unknown_not_scored',cases.loc[cases.experimental_label.isna(),'correct_when_determinate'].isna().all())
check('fixed_threshold_classes',(cases.predicted_label.eq(cases.score.gt(0).astype(int))).all())
contrasts=pd.read_csv(E/'pilot_paired_contrasts.csv')
check('contrast_pair_identity_flag_present','external_labels_identical_for_this_contrast' in contrasts)
check('contact_11_100_same_external_labels',contrasts.loc[contrasts.route.eq('contact') & contrasts.comparison.eq('AUROC11-AUROC100'),'external_labels_identical_for_this_contrast'].all())
check('oral_common_all_label_vectors_identical',contrasts.loc[contrasts.route.eq('oral'),'external_labels_identical_across_three_thresholds'].all())
for r in metrics.itertuples():
    d=cases[cases.route.eq(r.route)&cases['split'].eq(r.split)]
    if r.quality_screen=='explicit_adult_only':
        d=d[d.life_stage_evidence.eq('explicit_in_source')]
    if r.cohort=='common_three_thresholds':
        known=d.groupby('compound_id').experimental_label.apply(lambda x:x.notna().all())
        d=d[d.compound_id.isin(known.index[known])]
    d=d[d.threshold.eq(r.threshold)&d.experimental_label.notna()]
    check(f'cohort_n_{r.Index}',len(d)==r.n_compounds and d.compound_id.is_unique)
    if d.experimental_label.nunique()==2:
        check(f'scalar_auc_{r.Index}',abs(roc_auc_score(d.experimental_label,d.score)-r.AUROC)<1e-12)
        check(f'scalar_ap_{r.Index}',abs(average_precision_score(d.experimental_label,d.score)-r.average_precision)<1e-12)

report={'round':1,'status':'PASS_WITH_SCIENTIFIC_LIMITATIONS','checks_passed':len(checks),
        'checks':checks,'weighted_formula_sklearn_max_abs_error':float(max_error),
        'models_fitted_in_verifier':0,'note':'All pilot cohorts remain sparse; passing arithmetic checks is not adequate population validation.'}
(E/'round1_metric_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(f'First-round metric checks passed: {len(checks)}; maximum weighted AUROC error {max_error:.3g}.')
