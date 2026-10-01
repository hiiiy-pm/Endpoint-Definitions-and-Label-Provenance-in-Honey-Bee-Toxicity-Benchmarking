"""Evaluate frozen new-compound predictions after endpoint curation is frozen.

No model fitting or threshold optimisation. All small-sample metrics remain
exploratory; single-class discrimination metrics are explicitly undefined.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import beta
from sklearn.metrics import roc_auc_score, average_precision_score

OUT=Path(__file__).resolve().parent
DATA=OUT/'data_curation'
RESULTS=OUT/'evaluation'
RESULTS.mkdir(exist_ok=True)
THRESHOLDS=[1,11,100]
B=5000
SEED=20260923

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def truth(value,qualifier,t):
    v=float(value); q=str(qualifier).strip()
    if q in ('','=','nan'):
        return int(v<=t)
    if q=='>' and v>=t:
        return 0
    if q=='>=' and v>t:
        return 0
    if q in ('<','<=') and v<=t:
        return 1
    return np.nan

def exact_binomial(k,n):
    if n==0:
        return (np.nan,np.nan)
    return (0. if k==0 else float(beta.ppf(.025,k,n-k+1)),
            1. if k==n else float(beta.ppf(.975,k+1,n-k)))

def weighted_auc_many(y,s,w):
    pos=np.flatnonzero(y==1); neg=np.flatnonzero(y==0)
    out=np.full(len(w),np.nan)
    if not len(pos) or not len(neg):
        return out
    pw=w[:,pos]; nw=w[:,neg]
    comparison=(s[pos,None]>s[None,neg]).astype(float)+.5*(s[pos,None]==s[None,neg])
    numerator=np.einsum('bi,ij,bj->b',pw,comparison,nw)
    den=pw.sum(axis=1)*nw.sum(axis=1)
    np.divide(numerator,den,out=out,where=den>0)
    return out

def ci(values):
    finite=values[np.isfinite(values)]
    return (float(np.quantile(finite,.025)),float(np.quantile(finite,.975)),len(finite)) if len(finite) else (np.nan,np.nan,0)

def main():
    started=datetime.now(timezone.utc).isoformat()
    fm=json.loads((DATA/'endpoint_freeze_manifest.json').read_text(encoding='utf-8'))
    for name,digest in fm['files'].items():
        assert sha(DATA/name)==digest,name
    prediction_paths=[OUT/'predictions'/batch/'new_compound_predictions.csv' for batch in ['initial8','added2']]
    for path in prediction_paths:
        meta=json.loads((path.parent/'prediction_run_manifest.json').read_text(encoding='utf-8'))
        assert meta['completed'] and sha(path)==meta['prediction_sha256']
    predictions=pd.concat([pd.read_csv(p) for p in prediction_paths],ignore_index=True)
    assert len(predictions)==90
    assert not predictions.duplicated(['compound_id','split','threshold_ug_bee']).any()
    assert predictions.passes_conservative_full_development_exclusion.all()
    assert not predictions.score_is_probability.any()
    assert ((predictions.score>0).astype(int)==predictions.predicted_label).all()
    endpoints=pd.read_csv(DATA/'frozen_scoreable_48h_endpoints.csv')
    assert not endpoints.duplicated(['compound_id','route']).any()
    assert endpoints.observation_h.eq(48).all()
    assert endpoints.endpoint_name.eq('LD50').all()
    assert endpoints.species.eq('Apis mellifera').all()
    assert endpoints.compound_id.isin(predictions.compound_id).all()
    assert not endpoints[['overlap_inchikey','overlap_connectivity','overlap_parent_connectivity']].any().any()
    for t in THRESHOLDS:
        expected=np.array([truth(v,q,t) for v,q in zip(endpoints.LD50_value,endpoints.qualifier)])
        stored=pd.to_numeric(endpoints[f'label_le_{t}'],errors='coerce').to_numpy()
        assert np.allclose(expected,stored,equal_nan=True),f'interval truth {t}'
        endpoints[f'y_{t}']=expected
    joined=[]
    for _,e in endpoints.iterrows():
        pred=predictions[predictions.compound_id.eq(e.compound_id)]
        assert len(pred)==9
        for _,r in pred.iterrows():
            t=int(r.threshold_ug_bee)
            joined.append({'endpoint_id':e.endpoint_id,'compound_id':e.compound_id,'route':e.route,
                'observation_h':e.observation_h,'life_stage_evidence':e.life_stage_evidence,
                'qualifier':e.qualifier,'LD50_value':e.LD50_value,'original_unit':e.original_unit,
                'study_id':e.study_id,'source_url':e.source_url,'source_location':e.source_location,
                'split':r['split'],'threshold':t,'score':r.score,'predicted_label':r.predicted_label,
                'experimental_label':e[f'y_{t}'],'nearest_training_tanimoto':r.nearest_training_tanimoto,
                'model_artifact_sha256':r.model_artifact_sha256})
    joined=pd.DataFrame(joined)
    joined['correct_when_determinate']=np.where(joined.experimental_label.isna(),np.nan,
        (joined.predicted_label==joined.experimental_label).astype(int))
    joined.to_csv(RESULTS/'new_compound_case_predictions.csv',index=False,encoding='utf-8-sig')
    rows=[]; contrasts=[]; coverage=[]; boot_payload={}
    for ri,route in enumerate(['contact','oral']):
        route_data=endpoints[endpoints.route.eq(route)]
        for qi,quality in enumerate(['regulatory_context','explicit_adult_only']):
            subset=route_data if quality=='regulatory_context' else route_data[route_data.life_stage_evidence.eq('explicit_in_source')]
            for ci_i,cohort in enumerate(['threshold_specific','common_three_thresholds']):
                common=subset[subset[[f'y_{t}' for t in THRESHOLDS]].notna().all(axis=1)]
                for si,split in enumerate(['Random','MaxMin','Time']):
                    auc_boot={}; auc_point={}; sizes={}
                    if cohort=='common_three_thresholds':
                        fixed=common.sort_values('compound_id')
                        rng=np.random.default_rng(SEED+1000*ri+100*qi+10*ci_i+si)
                        weights=rng.multinomial(len(fixed),np.full(len(fixed),1/len(fixed)),size=B) if len(fixed) else np.zeros((B,0))
                    for ti,t in enumerate(THRESHOLDS):
                        selected=common if cohort=='common_three_thresholds' else subset[subset[f'y_{t}'].notna()]
                        selected=selected.sort_values('compound_id')
                        pred=predictions[predictions['split'].eq(split)&predictions.threshold_ug_bee.eq(t)].set_index('compound_id')
                        ids=selected.compound_id.tolist()
                        y=selected[f'y_{t}'].to_numpy(dtype=int)
                        s=pred.loc[ids,'score'].to_numpy()
                        binary=pred.loc[ids,'predicted_label'].to_numpy(dtype=int)
                        n=len(y); np_=int(y.sum()); nn=n-np_
                        tp=int(((binary==1)&(y==1)).sum()); tn=int(((binary==0)&(y==0)).sum())
                        fp=nn-tn; fn=np_-tp
                        point=float(roc_auc_score(y,s)) if np_ and nn else np.nan
                        ap=float(average_precision_score(y,s)) if np_ and nn else np.nan
                        if cohort=='threshold_specific':
                            rng=np.random.default_rng(SEED+1000*ri+100*qi+10*ci_i+si+10000*ti)
                            weights=rng.multinomial(n,np.full(n,1/n),size=B) if n else np.zeros((B,0))
                        dist=weighted_auc_many(y,s,weights)
                        lo,hi,valid=ci(dist)
                        auc_point[t]=point; auc_boot[t]=dist; sizes[t]=(np_,nn)
                        sens_ci=exact_binomial(tp,np_); spec_ci=exact_binomial(tn,nn)
                        rows.append({'route':route,'quality_screen':quality,'cohort':cohort,'split':split,
                            'threshold':t,'n_compounds':n,'n_positive':np_,'n_negative':nn,
                            'excluded_unknown_in_source_pool':int(subset[f'y_{t}'].isna().sum()),
                            'TP':tp,'TN':tn,'FP':fp,'FN':fn,
                            'observed_agreement':(tp+tn)/n if n else np.nan,
                            'sensitivity':tp/np_ if np_ else np.nan,'sensitivity_exact95_low':sens_ci[0],'sensitivity_exact95_high':sens_ci[1],
                            'specificity':tn/nn if nn else np.nan,'specificity_exact95_low':spec_ci[0],'specificity_exact95_high':spec_ci[1],
                            'AUROC':point,'average_precision':ap,'AUROC_exploratory95_low':lo,'AUROC_exploratory95_high':hi,
                            'bootstrap_valid':valid,'bootstrap_degenerate':B-valid,'bootstrap_requested':B,
                            'exploratory_sparse':min(np_,nn)<5,
                            'single_class':not (np_ and nn),'sampling_note':'Evidence-limited selected new-chemical pilot; not population prevalence'})
                        if cohort=='common_three_thresholds':
                            boot_payload[f'{route}__{quality}__{split}__auc{t}']=dist
                    if cohort=='common_three_thresholds':
                        labels_identical=len(common)>0 and common[['y_1','y_11','y_100']].nunique(axis=1).eq(1).all()
                        for high in [11,1]:
                            delta=auc_boot[high]-auc_boot[100]
                            lo,hi,valid=ci(delta)
                            point=auc_point[high]-auc_point[100]
                            contrasts.append({'route':route,'quality_screen':quality,'split':split,
                                'n_common_compounds':len(common),'comparison':f'AUROC{high}-AUROC100',
                                'estimate':point,'exploratory95_low':lo,'exploratory95_high':hi,
                                'bootstrap_valid':valid,'bootstrap_degenerate':B-valid,
                                'external_labels_identical_across_three_thresholds':labels_identical,
                                'external_labels_identical_for_this_contrast':bool(len(common)>0 and common[f'y_{high}'].eq(common['y_100']).all()),
                                'exploratory_sparse':min(*sizes[high],*sizes[100])<5,
                                'interpretation':'Different frozen-model ranking on a pilot; not a replication of label-threshold effects if intermediate tiers are absent'})
                            boot_payload[f'{route}__{quality}__{split}__delta{high}-100']=delta
            # Four-tier assessment without imputation of unknown boundaries.
            for row in subset.itertuples():
                vals=(getattr(row,'y_1'),getattr(row,'y_11'),getattr(row,'y_100'))
                tier='unresolved' if any(pd.isna(x) for x in vals) else str(sum(vals))
                coverage.append({'route':route,'quality_screen':quality,'compound_id':row.compound_id,'tier_0_to_3':tier})
    pd.DataFrame(rows).to_csv(RESULTS/'pilot_performance.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(contrasts).to_csv(RESULTS/'pilot_paired_contrasts.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(coverage).to_csv(RESULTS/'pilot_four_tier_members.csv',index=False,encoding='utf-8-sig')
    np.savez_compressed(RESULTS/'pilot_bootstrap_distributions.npz',**boot_payload)
    nested=[]
    for (compound,split),g in predictions.groupby(['compound_id','split']):
        d=g.set_index('threshold_ug_bee').predicted_label
        nested.append({'compound_id':compound,'split':split,'predicted_le_1':int(d[1]),
            'predicted_le_11':int(d[11]),'predicted_le_100':int(d[100]),
            'nested_class_order_violated':not(d[1]<=d[11]<=d[100])})
    pd.DataFrame(nested).to_csv(RESULTS/'predicted_class_nesting_diagnostic.csv',index=False,encoding='utf-8-sig')
    verify={'evaluation_started_utc':started,'completed_utc':datetime.now(timezone.utc).isoformat(),
        'endpoint_freeze_manifest_sha256':sha(DATA/'endpoint_freeze_manifest.json'),
        'prediction_sha256':{str(p.relative_to(OUT)):sha(p) for p in prediction_paths},
        'n_predicted_structures':int(predictions.compound_id.nunique()),'n_model_predictions':len(predictions),
        'n_scoreable_endpoints':len(endpoints),'n_case_model_rows':len(joined),
        'external_structures_used_in_fit':0,'fit_calls_in_evaluation':0,
        'all_endpoint_qualifier_labels_recalculated_and_equal':True,
        'no_single_class_AUROC_or_AP':bool(pd.DataFrame(rows).loc[lambda d:d.single_class,['AUROC','average_precision']].isna().all().all()),
        'all_cases_retained_without_score_selection':True,'script_sha256':sha(__file__),
        'bootstrap_seed':SEED,'bootstrap_repetitions':B}
    (RESULTS/'round1_evaluation_checks.json').write_text(json.dumps(verify,indent=2),encoding='utf-8')
    print('New-compound pilot evaluated. Candidate counts and all outcomes retained. No fits performed.')
    print(pd.DataFrame(rows).query("quality_screen=='regulatory_context' and cohort=='threshold_specific'")[['route','split','threshold','n_compounds','n_positive','n_negative','TP','TN','FP','FN','AUROC','average_precision']].to_string(index=False))

if __name__=='__main__':
    main()
