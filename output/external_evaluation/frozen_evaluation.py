"""Post hoc external-source reassessment of SAVED predictions; never fits a model.

Run with Python, pandas and NumPy. Input files are read-only and fingerprinted.
Same-compound paired bootstrap weights are shared across endpoints and label sets.
This is NOT new-compound independent external validation: these are original test
compounds with external-source endpoint labels and unresolved study provenance.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import importlib.util
import json
import re
import sys

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
AUDIT = OUT.parent / 'external_data_audit'
FEAS = OUT.parent / 'frozen_external_feasibility'
THRESHOLDS = (100, 11, 1)
SPLITS = ('random', 'maxmin', 'time')
REPS = ('ECFP', 'Avalon', 'MACCS', 'WL-HI', 'Descriptors12')
SOURCE_NAMES = ('OFT_single_structure_contact', 'OFT_48h_contact', 'PLOS_resolved_contact')
B = 5000
SEED = 20260922
CHECKS = {}


def check(name, condition):
    CHECKS[name] = bool(condition)
    if not condition:
        raise AssertionError(name)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(df, name):
    df.to_csv(OUT / name, index=False, encoding='utf-8-sig', float_format='%.17g')


def norm(x):
    return re.sub('[^a-z0-9]', '', str(x).lower())


def norm_name(x):
    return norm(re.sub(r'\s*\(Ref:.*', '', str(x), flags=re.I))


def text_value(x):
    return '' if pd.isna(x) else str(x)


def finite(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else None
    except (ValueError, TypeError):
        return None


def interval_label(lo, lower_qualifier, hi, upper_qualifier, threshold):
    """Match the prior documented audit's conservative interval convention."""
    lo, hi = finite(lo), finite(hi)
    lq, hq = text_value(lower_qualifier), text_value(upper_qualifier)
    if lo is not None and hi is None and lq in ('', '='):
        return int(lo <= threshold)
    if lo is not None and lq in ('<', '<=') and lo <= threshold:
        return 1
    if hi is not None and hq in ('<', '<=', '') and hi <= threshold:
        return 1
    if lo is not None and lq in ('>', '>=') and (lo > threshold or (lo == threshold and lq == '>')):
        return 0
    return -1


def identity_overlap(ext, train):
    """Exclude exact structures, connectivity, parent, CAS or normalized name."""
    matched = pd.Series(False, index=ext.index)
    for field in ('inchikey', 'connectivity', 'parent_connectivity', 'CAS', 'name'):
        transform = norm if field == 'CAS' else norm_name if field == 'name' else str
        values = {transform(x) for x in train[field].fillna('') if str(x)} - {''}
        matched |= ext[field].fillna('').map(transform).isin(values)
    return matched


def weighted_metrics(y, scores, weights):
    """Tie-aware weighted AUROC and step average precision, vectorized over B.

    AP = sum over distinct decreasing score groups of delta(recall)*precision.
    Both outputs are NaN when the resample contains only one outcome class.
    """
    y = np.asarray(y, dtype=int)
    scores = np.asarray(scores, dtype=float)
    weights = np.atleast_2d(np.asarray(weights, dtype=float))
    if len(y) == 0:
        return np.full(len(weights), np.nan), np.full(len(weights), np.nan)
    order = np.argsort(scores, kind='stable')
    sorted_score, sorted_y = scores[order], y[order]
    first = np.r_[0, np.flatnonzero(np.diff(sorted_score) != 0) + 1]
    w = weights[:, order]
    positive = np.add.reduceat(w * sorted_y, first, axis=1)
    negative = np.add.reduceat(w * (1 - sorted_y), first, axis=1)
    pos, neg = positive.sum(axis=1), negative.sum(axis=1)
    valid = (pos > 0) & (neg > 0)
    numerator = (positive * (np.cumsum(negative, axis=1) - .5 * negative)).sum(axis=1)
    auc = np.full(len(weights), np.nan)
    auc[valid] = numerator[valid] / (pos[valid] * neg[valid])
    p_desc, n_desc = positive[:, ::-1], negative[:, ::-1]
    cum_p = np.cumsum(p_desc, axis=1)
    cum_total = np.cumsum(p_desc + n_desc, axis=1)
    precision = np.divide(cum_p, cum_total, out=np.zeros_like(cum_p), where=cum_total > 0)
    ap = np.full(len(weights), np.nan)
    ap[valid] = (precision * p_desc).sum(axis=1)[valid] / pos[valid]
    return auc, ap


def brute_metrics(y, score):
    """Independent small-array references, without sorting/group reduction."""
    y, score = np.asarray(y), np.asarray(score)
    if len(np.unique(y)) != 2:
        return np.nan, np.nan
    positive, negative = score[y == 1], score[y == 0]
    auc = np.mean([(p > n) + .5 * (p == n) for p in positive for n in negative])
    ap = np.mean([y[score >= p].mean() for p in positive])
    return auc, ap


def metric_checks(predictions, official_metrics):
    examples = [([0, 1], [0, 1], 1., 1.), ([0, 1], [1, 0], 0., .5),
                ([0, 1], [1, 1], .5, .5), ([0, 0, 1, 1], [.1, .4, .35, .8], .75, 5 / 6)]
    for i, (y, s, expected_auc, expected_ap) in enumerate(examples):
        a, p = weighted_metrics(y, s, np.ones((1, len(y))))
        check(f'metric_known_example_{i}', np.allclose([a[0], p[0]], [expected_auc, expected_ap]))
    for positive_label in (0, 1):
        a, p = weighted_metrics([positive_label] * 4, [0, 1, 2, 3], np.ones((3, 4)))
        check(f'single_class_{positive_label}_both_metrics_NA', np.isnan(a).all() and np.isnan(p).all())
    a, p = weighted_metrics([0, 1], [0, 1], [[2, 0], [0, 2], [1, 1]])
    check('single_class_bootstrap_draws_NA', np.isnan(a[:2]).all() and np.isnan(p[:2]).all() and a[2] == 1 and p[2] == 1)
    rng = np.random.default_rng(7721)
    for i in range(100):
        y = np.r_[0, 1, rng.integers(0, 2, 6)]
        s = rng.integers(0, 4, len(y))  # Deliberately many tied scores.
        w = rng.multinomial(len(y), np.full(len(y), 1 / len(y)))
        a, p = weighted_metrics(y, s, w)
        expected = brute_metrics(np.repeat(y, w), np.repeat(s, w))
        check(f'metric_tie_bootstrap_enumeration_{i}', np.allclose([a[0], p[0]], expected, equal_nan=True))
    errors = []
    sklearn_available = importlib.util.find_spec('sklearn') is not None
    if sklearn_available:
        from sklearn.metrics import roc_auc_score, average_precision_score
    for (split, rep, threshold), g in predictions.groupby(['Split', 'Representation', 'Threshold']):
        a, p = weighted_metrics(g.Y, g.Score, np.ones((1, len(g))))
        row = official_metrics.loc[(official_metrics.Split == split) &
                                  (official_metrics.Representation == rep) &
                                  (official_metrics.Threshold == threshold)].iloc[0]
        errors.extend([abs(a[0] - row.AUROC), abs(p[0] - row.AUPRC)])
        if sklearn_available:
            check(f'sklearn_{split}_{rep}_{threshold}', np.allclose(
                [a[0], p[0]], [roc_auc_score(g.Y, g.Score), average_precision_score(g.Y, g.Score)], atol=1e-12))
    check('all_45_original_auc_and_ap_reproduced', max(errors) < 1e-12)
    return {'max_absolute_error_vs_45_original_metrics': float(max(errors)),
            'sklearn_independent_check_available': sklearn_available,
            'independent_small_tied_weighted_examples': 100, 'known_examples': 4}


def ci_fields(values):
    if values is None:
        return {'ci_low': np.nan, 'ci_high': np.nan, 'bootstrap_requested': 0,
                'bootstrap_valid': 0, 'bootstrap_degenerate': 0}
    finite_values = values[np.isfinite(values)]
    low, high = np.percentile(finite_values, [2.5, 97.5]) if len(finite_values) else (np.nan, np.nan)
    return {'ci_low': low, 'ci_high': high, 'bootstrap_requested': B,
            'bootstrap_valid': len(finite_values), 'bootstrap_degenerate': B - len(finite_values)}


def sparse_note(sparse, defined=True):
    if not defined:
        return 'Undefined: at least one endpoint/label set contains only one class.'
    if sparse:
        return ('Exploratory sparse classes (<5 in at least one class). Ordinary bootstrap excludes '
                'single-class draws; conditional percentile intervals may be misleading, especially '
                'with one positive. A narrow interval does not establish replication.')
    return 'Exploratory post hoc percentile interval; study provenance not demonstrated independent.'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    protocol = {
        'analysis': 'Frozen prediction external-source reassessment, designed after viewing labels',
        'not_independent_new_compound_validation': True,
        'fitting_calibration_tuning_or_score_transformation': False,
        'primary_representation_designated': 'ECFP',
        'other_representations': 'All four reported descriptively; none selected by performance.',
        'sources': {
            SOURCE_NAMES[0]: 'Audited OFT single_structure_screen AND route Contact',
            SOURCE_NAMES[1]: 'Audited OFT 48h_no_material_flag_screen AND route Contact; missing test material text is not proof of pure active ingredient',
            SOURCE_NAMES[2]: 'Audited PLOS resolved full InChIKey AND route Contact'},
        'cohorts': {
            'all_external_contact': 'Original mixed-route model reassessed against external contact endpoint',
            'original_contact_only': 'Fixed subset also labelled Contact in original ApisTox; other original study conditions can still differ',
            'quality_screen': 'PLOS only: no InChI/SMILES conflict in any selected record AND every selected record has an MRID; the implemented eligibility rule uses no prediction scores',
            'original_contact_quality_screen': 'PLOS only: intersection of original Contact and quality_screen'},
        'identity_exclusion': 'All members matched to model training set by full key, connectivity, standardized parent connectivity, CAS or normalized name removed; exclude entire external key if any row matches.',
        'identity_join': 'One-to-one exact full InChIKey to audited ApisTox, one-to-one SMILES to original processed data, verified CID and global Index and official split membership.',
        'external_labels': 'Audited censored intervals; every selected record for a key must have the same known label at each threshold; common members across all three thresholds.',
        'original_labels': 'Processed original Tier >=1, >=2, >=3 for 100,11,1 respectively, each checked against stored Y.',
        'contrasts': ['AUROC11 - AUROC100', 'AUROC1 - AUROC100', 'external minus original for each metric and each contrast'],
        'bootstrap': {'unit': 'compound', 'method': 'ordinary multinomial paired', 'B': B, 'seed': SEED,
                      'same_weights': 'Across all thresholds and original/external labels within each source/split/cohort',
                      'CI': 'Unadjusted percentile 95%; report valid and degenerate draws',
                      'single_class': 'NA; exclude undefined draws separately for each metric or contrast',
                      'sparse_flag': 'Any required positive or negative class count <5; no confirmatory interpretation'},
        'aggregation': 'Do not pool sources, splits, or overlapping cohorts as independent samples.',
        'quality_limitations': 'External regulatory data may share underlying studies with ApisTox; OFT broad duration/unit/material heterogeneity remains; original toxicity endpoint conditions and routes may differ.'
    }
    (OUT / 'frozen_protocol.json').write_text(json.dumps(protocol, ensure_ascii=False, indent=2), encoding='utf-8')
    inputs = [AUDIT / 'apistox_structures.csv', AUDIT / 'oft_acute_endpoints_aligned.csv',
              AUDIT / 'plos_adult_acute_endpoints_aligned.csv', FEAS / 'frozen_split_external_members.csv',
              ROOT / 'results/primary/02_test_predictions.csv', ROOT / 'results/primary/01_representation_performance.csv',
              ROOT / 'results/cache/processed_data.csv']
    inputs += [ROOT / f'data/official_splits/{s}_{part}.csv' for s in SPLITS for part in ('train', 'test')]
    inputs.append(FEAS / 'REBUILD_MANIFEST.json')
    initial_hashes = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    member_manifest = json.loads((FEAS / 'REBUILD_MANIFEST.json').read_text(encoding='utf-8'))
    check('membership_rebuild_manifest_passed', member_manifest['status'] == 'passed')
    check('membership_generator_hash', sha(ROOT / member_manifest['generator']) == member_manifest['generator_sha256'])
    for relative, expected_hash in member_manifest['input_sha256'].items():
        check(f'membership_input_hash_{relative}', sha(ROOT / relative) == expected_hash)
    for relative, expected_hash in member_manifest['output_sha256'].items():
        check(f'membership_output_hash_{relative}', sha(ROOT / relative) == expected_hash)
    a = pd.read_csv(inputs[0]).fillna('')
    oft = pd.read_csv(inputs[1]).fillna('')
    plos = pd.read_csv(inputs[2]).fillna('')
    prior_members = pd.read_csv(inputs[3])
    pred = pd.read_csv(inputs[4])
    official_metrics = pd.read_csv(inputs[5])
    processed = pd.read_csv(inputs[6])
    check('1035_unique_original_structures', len(a) == 1035 and a.inchikey.nunique() == 1035 and a.SMILES.nunique() == 1035)
    check('processed_structure_alignment', len(processed) == 1035 and processed.SMILES.nunique() == 1035)
    check('9315_predictions_45_cells', len(pred) == 9315 and pred.groupby(['Split', 'Representation', 'Threshold']).size().eq(207).all())
    check('one_score_per_original_index_cell', not pred.duplicated(['Split', 'Representation', 'Threshold', 'Index']).any())
    check('no_nonfinite_scores', np.isfinite(pred.Score).all())
    processed = processed.reset_index(names='Index').merge(a[['SMILES', 'inchikey', 'tier']], on='SMILES', validate='one_to_one')
    check('audit_original_tier_match', processed.Tier.eq(processed.tier).all())
    for t, cutoff in ((100, 1), (11, 2), (1, 3)):
        check(f'tier_to_original_Y{t}', processed[f'Y{t}'].eq(processed.Tier.ge(cutoff).astype(int)).all())
    lookup = processed.set_index('Index')
    check('prediction_CID_matches_original_index', pred.CID.eq(pred.Index.map(lookup.CID)).all())
    for t in THRESHOLDS:
        sub = pred[pred.Threshold.eq(t)]
        check(f'prediction_original_Y_{t}', sub.Y.eq(sub.Index.map(lookup[f'Y{t}'])).all())
    metric_verification = metric_checks(pred, official_metrics)
    for source, frame in [('OFT', oft), ('PLOS', plos)]:
        for t in THRESHOLDS:
            rebuilt = [interval_label(r.get('lower_ug_bee', ''), r.get('lower_qualifier', ''),
                                      r.get('upper_ug_bee', ''), r.get('upper_qualifier', ''), t)
                       for r in frame.to_dict('records')]
            check(f'{source}_audited_interval_labels_{t}', np.array_equal(rebuilt, frame[f'y_{t}'].to_numpy()))
    sets = {
        SOURCE_NAMES[0]: oft[oft.single_structure_screen & oft.route.eq('Contact')],
        SOURCE_NAMES[1]: oft[oft['48h_no_material_flag_screen'] & oft.route.eq('Contact')],
        SOURCE_NAMES[2]: plos[plos.route.eq('Contact') & plos.inchikey.ne('')],
    }
    members_rows, endpoint_rows, split_checks = [], [], []
    for split in SPLITS:
        train_file = pd.read_csv(ROOT / f'data/official_splits/{split}_train.csv')
        test_file = pd.read_csv(ROOT / f'data/official_splits/{split}_test.csv')
        train = a[a.SMILES.isin(train_file.SMILES)]
        test = a[a.SMILES.isin(test_file.SMILES)]
        check(f'{split}_official_sizes', len(train) == 828 and len(test) == 207)
        check(f'{split}_disjoint_official_SMILES', not set(train.SMILES) & set(test.SMILES))
        expected = processed[processed.SMILES.isin(test_file.SMILES)]
        s_pred = pred[pred.Split.str.lower().eq(split)]
        for (rep, t), g in s_pred.groupby(['Representation', 'Threshold']):
            check(f'{split}_{rep}_{t}_official_test_indexes', set(g.Index) == set(expected.Index))
            expected_order = test_file.reset_index(names='TestOrder').merge(processed[['SMILES', 'Index']], on='SMILES', validate='one_to_one')
            check(f'{split}_{rep}_{t}_official_test_order', g.sort_values('TestOrder').Index.tolist() == expected_order.Index.tolist())
        for source, external in sets.items():
            overlap = identity_overlap(external, train)
            excluded_keys = set(external.loc[overlap, 'inchikey'])
            eligible = external.loc[~overlap & ~external.inchikey.isin(excluded_keys)]
            check(f'{split}_{source}_no_training_identity_match', not identity_overlap(eligible, train).any())
            rebuilt = []
            for key, g in eligible.groupby('inchikey', sort=True):
                row = {'source': source, 'split': split, 'inchikey': key}
                for t in THRESHOLDS:
                    values = set(g[f'y_{t}'].astype(int))
                    row[f'y_{t}'] = next(iter(values)) if len(values) == 1 and -1 not in values else -1
                row['exact_match_to_original_test'] = key in set(test.inchikey)
                rebuilt.append(row)
            rebuilt = pd.DataFrame(rebuilt)
            prior = prior_members[prior_members.source.eq(source) & prior_members.split.eq(split)]
            cols = ['inchikey', 'y_1', 'y_11', 'y_100', 'exact_match_to_original_test']
            check(f'{split}_{source}_prior_members_independently_rebuilt',
                  rebuilt[cols].sort_values('inchikey').reset_index(drop=True).equals(prior[cols].sort_values('inchikey').reset_index(drop=True)))
            fixed = rebuilt[rebuilt.exact_match_to_original_test & rebuilt[['y_1', 'y_11', 'y_100']].ge(0).all(axis=1)].copy()
            fixed = fixed.merge(processed[['Index', 'CID', 'inchikey', 'name', 'SMILES', 'toxicity_type', 'Tier', 'Y100', 'Y11', 'Y1']], on='inchikey', validate='one_to_one')
            check(f'{split}_{source}_nested_external_labels', (fixed.y_1 <= fixed.y_11).all() and (fixed.y_11 <= fixed.y_100).all())
            for row in fixed.to_dict('records'):
                g = eligible[eligible.inchikey.eq(row['inchikey'])]
                row['n_external_records'] = len(g)
                row['n_exact_records'] = int(g.value_kind.eq('exact').sum())
                row['external_route'] = 'Contact'
                row['original_route'] = row.pop('toxicity_type')
                for col in ('duration_hours', 'unit_mapping', 'dose_basis', 'MRID', 'source_dois', 'source_years', 'MRID_mapping_status', 'mapping_method'):
                    row[f'external_{col}'] = '|'.join(sorted({str(v) for v in g[col] if str(v)})) if col in g else 'not_in_table'
                row['external_material_review_flag_any'] = bool(g.material_review_needed.any()) if 'material_review_needed' in g else 'not_in_table'
                row['external_inchi_smiles_disagree_any'] = bool(g.inchi_smiles_disagree.any())
                row['external_MRID_all_present'] = bool(g.MRID.map(lambda x: bool(re.search(r'\d{6,}', str(x)))).all()) if 'MRID' in g else False
                row['passes_plos_quality_screen'] = bool(not row['external_inchi_smiles_disagree_any'] and row['external_MRID_all_present']) if source.startswith('PLOS') else False
                row['minimum_lower_ug_bee'] = pd.to_numeric(g.lower_ug_bee, errors='coerce').min()
                row['record_values_with_qualifiers'] = '|'.join(sorted({text_value(r['lower_qualifier']) + str(r['lower_ug_bee']) for r in g.to_dict('records')}))
                members_rows.append(row)
                for record in g.to_dict('records'):
                    record.update(analysis_source=source, analysis_split=split, original_Index=row['Index'], original_CID=row['CID'])
                    endpoint_rows.append(record)
            split_checks.append({'source': source, 'split': split, 'external_records_before_train_filter': len(external),
                                 'train_identity_matching_records': int(overlap.sum()), 'remaining_eligible_records': len(eligible),
                                 'remaining_eligible_keys': eligible.inchikey.nunique(), 'fixed_three_endpoint_scored_compounds': len(fixed),
                                 'fixed_original_contact_compounds': int(fixed.toxicity_type.eq('Contact').sum())})
    members = pd.DataFrame(members_rows)
    save(members, 'frozen_compound_members.csv')
    save(pd.DataFrame(endpoint_rows), 'frozen_underlying_endpoint_rows.csv')
    save(pd.DataFrame(split_checks), 'frozen_cohort_flow.csv')
    rows, contrast_rows, change_rows, score_rows = [], [], [], []
    distributions = {}
    for source_i, source in enumerate(SOURCE_NAMES):
        for split_i, split in enumerate(SPLITS):
            base = members[members.source.eq(source) & members.split.eq(split)].sort_values('inchikey')
            cohorts = ['all_external_contact', 'original_contact_only']
            if source.startswith('PLOS'):
                cohorts += ['quality_screen', 'original_contact_quality_screen']
            for cohort_i, cohort in enumerate(cohorts):
                m = base
                if cohort in ('original_contact_only', 'original_contact_quality_screen'):
                    m = m[m.original_route.eq('Contact')]
                if cohort in ('quality_screen', 'original_contact_quality_screen'):
                    m = m[m.passes_plos_quality_screen]
                n = len(m)
                if n == 0:
                    continue
                rng = np.random.default_rng(SEED + source_i * 1000 + split_i * 100 + cohort_i)
                weights = rng.multinomial(n, np.full(n, 1 / n), size=B)
                context = {'source': source, 'split': split, 'cohort': cohort, 'n_compounds': n}
                point, boot, counts = {}, {}, {}
                for rep in REPS:
                    for t in THRESHOLDS:
                        subpred = pred[pred.Split.str.lower().eq(split) & pred.Representation.eq(rep) & pred.Threshold.eq(t)]
                        matched = m.merge(subpred[['Index', 'CID', 'Score', 'Y']], on=['Index', 'CID'], validate='one_to_one', how='left')
                        check(f'{source}_{split}_{cohort}_{rep}_{t}_all_scores_available', matched.Score.notna().all() and len(matched) == n)
                        score = matched.Score.to_numpy()
                        if cohort_i == 0:
                            score_rows.extend([{**context, 'Representation': rep, 'threshold': t,
                                                'inchikey': r['inchikey'], 'Index': r['Index'], 'CID': r['CID'],
                                                'Score': r['Score'], 'original_label': r[f'Y{t}'], 'external_label': r[f'y_{t}']}
                                               for r in matched.to_dict('records')])
                        for label_source, column in [('original', f'Y{t}'), ('external', f'y_{t}')]:
                            y = matched[column].astype(int).to_numpy()
                            pos = int(y.sum())
                            neg = n - pos
                            counts[(label_source, t)] = (pos, neg)
                            metric_vals = weighted_metrics(y, score, np.ones((1, n)))
                            metric_boot = weighted_metrics(y, score, weights) if rep == 'ECFP' else (None, None)
                            for metric, val, dist in zip(('AUROC', 'average_precision'), metric_vals, metric_boot):
                                key = (rep, label_source, t, metric)
                                point[key] = float(val[0])
                                boot[key] = dist
                                rows.append({**context, 'Representation': rep, 'threshold': t, 'labels': label_source,
                                             'metric': metric, 'estimate': val[0], 'n_pos': pos, 'n_neg': neg,
                                             'exploratory_sparse': min(pos, neg) < 5,
                                             'note': sparse_note(min(pos, neg) < 5, np.isfinite(val[0])), **ci_fields(dist)})
                                if dist is not None:
                                    distributions[f'{source}__{split}__{cohort}__{label_source}__{t}__{metric}'] = dist
                    for t in THRESHOLDS:
                        sparse = min(*counts[('original', t)], *counts[('external', t)]) < 5
                        for metric in ('AUROC', 'average_precision'):
                            delta = point[(rep, 'external', t, metric)] - point[(rep, 'original', t, metric)]
                            dist = boot[(rep, 'external', t, metric)] - boot[(rep, 'original', t, metric)] if rep == 'ECFP' else None
                            change_rows.append({**context, 'Representation': rep, 'threshold': t, 'metric': metric,
                                                'comparison': 'external_minus_original', 'estimate': delta,
                                                'original_n_pos': counts[('original', t)][0], 'external_n_pos': counts[('external', t)][0],
                                                'exploratory_sparse': sparse, 'note': sparse_note(sparse, np.isfinite(delta)), **ci_fields(dist)})
                    for high in (11, 1):
                        contrast_dist = {}
                        contrast_point = {}
                        for label_source in ('original', 'external'):
                            delta = point[(rep, label_source, high, 'AUROC')] - point[(rep, label_source, 100, 'AUROC')]
                            dist = boot[(rep, label_source, high, 'AUROC')] - boot[(rep, label_source, 100, 'AUROC')] if rep == 'ECFP' else None
                            contrast_dist[label_source] = dist
                            contrast_point[label_source] = delta
                            sparse = min(*counts[(label_source, high)], *counts[(label_source, 100)]) < 5
                            contrast_rows.append({**context, 'Representation': rep, 'labels': label_source,
                                                  'comparison': f'AUROC{high}-AUROC100', 'estimate': delta,
                                                  'n_pos_high_endpoint': counts[(label_source, high)][0],
                                                  'n_neg_high_endpoint': counts[(label_source, high)][1],
                                                  'n_pos_100': counts[(label_source, 100)][0], 'n_neg_100': counts[(label_source, 100)][1],
                                                  'exploratory_sparse': sparse, 'note': sparse_note(sparse, np.isfinite(delta)), **ci_fields(dist)})
                        shift = contrast_point['external'] - contrast_point['original']
                        dist = contrast_dist['external'] - contrast_dist['original'] if rep == 'ECFP' else None
                        sparse = min(*(v for lab in ('original', 'external') for t in (high, 100) for v in counts[(lab, t)])) < 5
                        contrast_rows.append({**context, 'Representation': rep, 'labels': 'external_minus_original',
                                              'comparison': f'AUROC{high}-AUROC100', 'estimate': shift,
                                              'exploratory_sparse': sparse, 'note': sparse_note(sparse, np.isfinite(shift)), **ci_fields(dist)})
    metrics = pd.DataFrame(rows)
    contrasts = pd.DataFrame(contrast_rows)
    score_table = pd.DataFrame(score_rows)
    for (source, split), table in score_table.groupby(['source', 'split']):
        comparison = table.merge(pred[pred.Split.str.lower().eq(split)],
                                 left_on=['Representation', 'threshold', 'Index', 'CID'],
                                 right_on=['Representation', 'Threshold', 'Index', 'CID'],
                                 validate='one_to_one', suffixes=('_export', '_saved'))
        check(f'{source}_{split}_exported_scores_exactly_frozen', len(comparison) == len(table) and comparison.Score_export.eq(comparison.Score_saved).all())
        check(f'{source}_{split}_exported_original_labels_exactly_frozen', comparison.original_label.eq(comparison.Y).all())
    ecfp_metric = metrics[metrics.Representation.eq('ECFP')]
    check('metric_bootstrap_valid_plus_degenerate_equals_B', (ecfp_metric.bootstrap_valid + ecfp_metric.bootstrap_degenerate).eq(B).all())
    check('point_metric_missing_iff_single_class', metrics.estimate.isna().eq(metrics.n_pos.eq(0) | metrics.n_neg.eq(0)).all())
    check('single_class_no_valid_bootstrap', ecfp_metric.loc[ecfp_metric.estimate.isna(), 'bootstrap_valid'].eq(0).all())
    check('metric_estimates_in_unit_interval', metrics.estimate.dropna().between(0, 1).all())
    check('metric_CI_in_unit_interval', metrics.ci_low.dropna().between(0, 1).all() and metrics.ci_high.dropna().between(0, 1).all())
    save(metrics, 'frozen_metrics_all_representations.csv')
    save(metrics[metrics.Representation.eq('ECFP')], 'frozen_metrics_ECFP.csv')
    save(contrasts, 'frozen_contrasts_all_representations.csv')
    save(contrasts[contrasts.Representation.eq('ECFP')], 'frozen_contrasts_ECFP.csv')
    save(pd.DataFrame(change_rows), 'frozen_external_minus_original_metrics.csv')
    save(score_table, 'frozen_scores_and_labels.csv')
    save(official_metrics, 'frozen_original_full_test_performance_context.csv')
    count_table = metrics[(metrics.Representation == 'ECFP') & (metrics.metric == 'AUROC')]
    save(count_table[['source', 'split', 'cohort', 'threshold', 'labels', 'n_compounds', 'n_pos', 'n_neg', 'exploratory_sparse']], 'frozen_all_cohort_counts.csv')
    np.savez_compressed(OUT / 'frozen_ECFP_paired_bootstrap_distributions.npz', **distributions)
    overlap_rows = []
    groups = list(members.groupby(['source', 'split']))
    for i, ((src1, spl1), g1) in enumerate(groups):
        for (src2, spl2), g2 in groups[i:]:
            overlap_rows.append({'source_1': src1, 'split_1': spl1, 'source_2': src2, 'split_2': spl2,
                                 'shared_compounds': len(set(g1.inchikey) & set(g2.inchikey))})
    save(pd.DataFrame(overlap_rows), 'frozen_cohort_overlap.csv')
    check('input_hashes_unchanged_after_analysis', all(sha(ROOT / p) == value for p, value in initial_hashes.items()))
    check('all_three_endpoint_counts_same_per_fixed_cohort', metrics.groupby(['source', 'split', 'cohort']).n_compounds.nunique().eq(1).all())
    check('zero_fitted_models', True)
    check('all_nonfinite_endpoint_estimates_marked', metrics[metrics.estimate.isna()].note.str.startswith('Undefined:').all())
    verification = {'analysis_complete': True, 'checks_passed': sum(CHECKS.values()), 'checks': CHECKS,
                    'metric_verification': metric_verification, 'input_sha256': initial_hashes,
                    'models_fitted': 0, 'predictions_recomputed': 0, 'bootstrap_repetitions': B,
                    'evaluation_script_sha256': sha(Path(__file__)),
                    'numpy_version': np.__version__, 'pandas_version': pd.__version__, 'python_version': sys.version,
                    'cohort_compound_counts': pd.DataFrame(split_checks).to_dict('records'),
                    'all_cohorts_reuse_original_test_compounds': True,
                    'unique_original_compounds_across_all_sources_and_splits': int(members.inchikey.nunique()),
                    'important_limitation': 'Post hoc external-label reassessment; data provenance may overlap and no new-compound independent external validation is claimed.'}
    verification['verification_round'] = 1
    (OUT / 'round1_verification.json').write_text(json.dumps(verification, ensure_ascii=False, indent=2), encoding='utf-8')
    (OUT / 'verification.json').write_text(json.dumps(verification, ensure_ascii=False, indent=2), encoding='utf-8')
    print('All input files unchanged; no models fitted. Checks:', len(CHECKS))
    main_auc = metrics[(metrics.Representation == 'ECFP') & (metrics.metric == 'AUROC') & (metrics.cohort == 'all_external_contact')]
    print(main_auc[['source', 'split', 'threshold', 'labels', 'n_compounds', 'n_pos', 'estimate', 'ci_low', 'ci_high', 'bootstrap_valid']].to_string(index=False))
    print('\nSAME ORIGINAL CONTACT SENSITIVITY:')
    sens = metrics[(metrics.Representation == 'ECFP') & (metrics.metric == 'AUROC') & (metrics.cohort == 'original_contact_only')]
    print(sens[['source', 'split', 'threshold', 'labels', 'n_compounds', 'n_pos', 'estimate', 'ci_low', 'ci_high']].to_string(index=False))


if __name__ == '__main__':
    main()
