"""Supplementary sparse descriptors with training-fold-only preprocessing.

Five validation folds select logistic C separately for each prespecified k.
Every fold fits its scaler and supervised stability ranking on that fold's
training rows. Elastic-net parameters use another five-fold CV within those
training rows. One fold's selector is cached for all C/k candidates. The final
selector uses the full outer training set; all fixed k values reach outer test.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import platform
import time
import warnings

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from _common import DESC_NAMES, RESULTS, ROOT, SEED, THRESHOLDS, build_or_load_cache

OUT = RESULTS / 'supplementary'
PROTOCOL = 'nested_training_fold_selector'
ALPHA_GRID = (1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2)
L1_GRID = (0.25, 0.5, 0.75)
C_GRID = (0.03, 0.1, 0.3, 1, 3, 10)
TOPK_GRID = (1, 2, 3, 5, 8, 12)
STABILITY_REPEATS = 500
N_FOLDS = 5
TABLES = ('01_elasticnet_stability_selection.csv',
          '02_topk_descriptor_performance.csv', '03_descriptor_consensus.csv')


def digest_array(values):
    return hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@dataclass
class Selector:
    scaler: StandardScaler
    order: np.ndarray
    frequency: np.ndarray
    median_coefficient: np.ndarray
    alpha: float
    l1_ratio: float
    tuning_auc: float


class FitAudit:
    """Check the row scope of each scaler and estimator before fitting."""

    def __init__(self):
        self.rows = []
        self.convergence_warnings = 0

    def fit(self, estimator, x, y, fit_ids, allowed_ids, context, kind):
        fit_ids, allowed_ids = np.asarray(fit_ids), np.asarray(allowed_ids)
        if len(x) != len(fit_ids) or (y is not None and len(y) != len(fit_ids)):
            raise AssertionError('Fit data and row identifiers differ in length.')
        forbidden_count = int(np.sum(~np.isin(fit_ids, allowed_ids)))
        if forbidden_count:
            raise AssertionError(f'Fit outside permitted rows: {context}/{kind}')
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always', ConvergenceWarning)
            warnings.filterwarnings('ignore', category=FutureWarning)
            estimator.fit(x, y) if y is not None else estimator.fit(x)
        convergence = sum(issubclass(w.category, ConvergenceWarning) for w in captured)
        self.convergence_warnings += convergence
        self.rows.append({'context': context, 'kind': kind, 'n_fit_rows': len(fit_ids),
                          'n_unique_fit_rows': len(np.unique(fit_ids)),
                          'n_allowed_rows': len(np.unique(allowed_ids)),
                          'n_forbidden_fit_rows': forbidden_count,
                          'ordered_fit_indices_sha256': digest_array(fit_ids.astype(np.int32)),
                          'allowed_indices_sha256': digest_array(np.sort(allowed_ids).astype(np.int32)),
                          'convergence_warnings': convergence})
        return estimator


def elasticnet_model(alpha, l1_ratio, random_state=SEED):
    return SGDClassifier(loss='log_loss', penalty='elasticnet', alpha=alpha,
                         l1_ratio=l1_ratio, class_weight='balanced', max_iter=1000,
                         tol=1e-4, random_state=random_state)


def logistic_model(c):
    return LogisticRegression(penalty='l2', solver='liblinear', C=c,
                              class_weight='balanced', max_iter=2000,
                              random_state=SEED)


def fit_selector(x_train, y_train, global_ids, sample_seed, audit, context,
                 repeats=STABILITY_REPEATS):
    """Fit a selector using only supplied training rows and their labels."""
    global_ids = np.asarray(global_ids)
    inner_cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
    scaled_folds = []
    for fold, (fit, val) in enumerate(inner_cv.split(x_train, y_train)):
        if set(global_ids[fit]) & set(global_ids[val]):
            raise AssertionError('Selector parameter CV folds overlap.')
        scaler = audit.fit(StandardScaler(), x_train[fit], None, global_ids[fit],
                           global_ids[fit], context, f'elasticnet_tuning_scaler_fold{fold}')
        scaled_folds.append((fit, val, scaler.transform(x_train[fit]),
                             scaler.transform(x_train[val])))
    best = (-np.inf, None, None)
    for alpha in ALPHA_GRID:
        for l1_ratio in L1_GRID:
            scores = []
            for fold, (fit, val, x_fit, x_val) in enumerate(scaled_folds):
                model = audit.fit(elasticnet_model(alpha, l1_ratio), x_fit, y_train[fit],
                                  global_ids[fit], global_ids[fit], context,
                                  f'elasticnet_tuning_model_fold{fold}')
                scores.append(roc_auc_score(y_train[val], model.decision_function(x_val)))
            mean_auc = float(np.mean(scores))
            if mean_auc > best[0]:  # Grid order resolves exact ties.
                best = (mean_auc, alpha, l1_ratio)
    scaler = audit.fit(StandardScaler(), x_train, None, global_ids, global_ids,
                       context, 'stability_scaler')
    standardized = scaler.transform(x_train)
    pos, neg = np.flatnonzero(y_train == 1), np.flatnonzero(y_train == 0)
    if min(len(pos), len(neg)) < N_FOLDS:
        raise AssertionError('Insufficient class counts for five-fold selector CV.')
    rng = np.random.default_rng(sample_seed)
    coefficients = np.zeros((repeats, x_train.shape[1]), dtype=float)
    for repeat in range(repeats):
        sample = np.concatenate([rng.choice(pos, len(pos), replace=True),
                                 rng.choice(neg, len(neg), replace=True)])
        rng.shuffle(sample)
        model = audit.fit(elasticnet_model(best[1], best[2], SEED + repeat),
                          standardized[sample], y_train[sample], global_ids[sample],
                          global_ids, context, 'stability_bootstrap')
        coefficients[repeat] = model.coef_[0]
    frequency = np.mean(np.abs(coefficients) > 1e-8, axis=0)
    median = np.median(coefficients, axis=0)
    order = np.lexsort((-np.abs(median), -frequency))
    return Selector(scaler, order, frequency, median, float(best[1]), float(best[2]), best[0])


def run_analysis(output_dir=OUT, audit_dir=None):
    started = time.perf_counter()
    output_dir = Path(output_dir)
    audit_dir = Path(audit_dir) if audit_dir else ROOT / 'output/descriptor_audit'
    output_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    cache = build_or_load_cache(force=False)
    study = cache['study']
    x = np.asarray(cache['descriptors'], dtype=float)
    audit = FitAudit()
    stability_rows, topk_rows, fold_rows, tuning_rows, predictions = [], [], [], [], []
    original_paths = [ROOT / 'data/raw/apistox.csv', RESULTS / 'cache/descriptors12.npy',
                      RESULTS / 'cache/processed_data.csv', ROOT / 'code/_common.py']
    original_paths += [ROOT / f'data/official_splits/{s}_{p}.csv'
                       for s in ('random', 'maxmin', 'time') for p in ('train', 'test')]
    original_hashes = {str(p.relative_to(ROOT)): file_hash(p) for p in original_paths}
    with threadpool_limits(limits=1):
        for split_i, (split, (train, test)) in enumerate(study.splits.items()):
            if len(train) != 828 or len(test) != 207 or set(train) & set(test):
                raise AssertionError('Unexpected original outer split.')
            for threshold_i, threshold in enumerate(THRESHOLDS):
                print(f'{split} threshold={threshold}: preparing five training-fold selectors', flush=True)
                y = study.y[threshold]
                x_train, y_train = x[train], y[train]
                base_seed = SEED + split_i * 1000 + threshold_i * 100
                validation_cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
                folds = []
                for fold, (fit, val) in enumerate(validation_cv.split(x_train, y_train)):
                    fit_global, val_global = train[fit], train[val]
                    if set(fit_global) & (set(val_global) | set(test)):
                        raise AssertionError('Fitting fold overlaps validation or outer test.')
                    context = f'{split}/threshold{threshold}/selection_fold{fold}'
                    selector = fit_selector(x_train[fit], y_train[fit], fit_global,
                                            base_seed + 10000 * (fold + 1), audit, context)
                    x_fit = selector.scaler.transform(x_train[fit])
                    x_val = selector.scaler.transform(x_train[val])
                    folds.append((fit, val, x_fit, x_val, selector))
                    fold_rows.append({'Split': split, 'Threshold': threshold, 'Fold': fold,
                                      'n_fit': len(fit), 'n_validation': len(val),
                                      'fit_indices': '|'.join(map(str, fit_global)),
                                      'validation_indices': '|'.join(map(str, val_global)),
                                      'outer_test_indices': '|'.join(map(str, test)),
                                      'BestAlpha': selector.alpha, 'BestL1Ratio': selector.l1_ratio,
                                      'SelectorCV_AUROC': selector.tuning_auc,
                                      'DescriptorOrder': ','.join(DESC_NAMES[j] for j in selector.order),
                                      'ScalerMean': '|'.join(map(repr, selector.scaler.mean_)),
                                      'ScalerScale': '|'.join(map(repr, selector.scaler.scale_))})
                # All C/k candidates share a fold's selector, never its labels.
                selected_c = {}
                for k in TOPK_GRID:
                    best = (-np.inf, None)
                    for c in C_GRID:
                        fold_scores = []
                        for fold, (fit, val, x_fit, x_val, selector) in enumerate(folds):
                            features = selector.order[:k]
                            model = audit.fit(logistic_model(c), x_fit[:, features], y_train[fit],
                                              train[fit], train[fit],
                                              f'{split}/threshold{threshold}/selection_fold{fold}',
                                              f'logistic_C_tuning_k{k}')
                            fold_scores.append(roc_auc_score(y_train[val], model.decision_function(x_val[:, features])))
                        mean_auc = float(np.mean(fold_scores))
                        tuning_rows.append({'Split': split, 'Threshold': threshold, 'TopK': k,
                                            'C': c, 'MeanInnerCV_AUROC': mean_auc,
                                            **{f'Fold{i}_AUROC': float(a) for i, a in enumerate(fold_scores)}})
                        if mean_auc > best[0]:
                            best = (mean_auc, c)
                    selected_c[k] = best
                final = fit_selector(x_train, y_train, train, base_seed, audit,
                                     f'{split}/threshold{threshold}/final_outer_training')
                x_final_train = final.scaler.transform(x_train)
                x_final_test = final.scaler.transform(x[test])
                for j, name in enumerate(DESC_NAMES):
                    stability_rows.append({'Split': split, 'Threshold': threshold, 'Descriptor': name,
                                           'SelectionFrequency': float(final.frequency[j]),
                                           'MedianCoefficient': float(final.median_coefficient[j]),
                                           'BestAlpha': final.alpha, 'BestL1Ratio': final.l1_ratio,
                                           'TrainCV_AUROC': final.tuning_auc, 'Protocol': PROTOCOL,
                                           'StabilityRepeats': STABILITY_REPEATS})
                for k in TOPK_GRID:
                    features = final.order[:k]
                    best_auc, best_c = selected_c[k]
                    model = audit.fit(logistic_model(best_c), x_final_train[:, features], y_train,
                                      train, train, f'{split}/threshold{threshold}/final_outer_training',
                                      f'logistic_final_k{k}')
                    scores = model.decision_function(x_final_test[:, features])
                    auc = float(roc_auc_score(y[test], scores))
                    topk_rows.append({'Split': split, 'Threshold': threshold, 'TopK': k,
                                      'TestAUROC': auc, 'SelectedDescriptors': ','.join(DESC_NAMES[j] for j in features),
                                      'BestC': best_c, 'InnerCV_AUROC': best_auc, 'Protocol': PROTOCOL,
                                      'InnerFolds': N_FOLDS, 'SelectorTuningFolds': N_FOLDS,
                                      'StabilityRepeats': STABILITY_REPEATS, 'n_train': len(train), 'n_test': len(test)})
                    predictions.extend({'Split': split, 'Threshold': threshold, 'TopK': k,
                                        'Index': int(index), 'CID': int(study.df.iloc[index].CID),
                                        'Y': int(y[index]), 'Score': float(score)}
                                       for index, score in zip(test, scores))
                print(f'{split} threshold={threshold}: completed all six k values', flush=True)
    stability = pd.DataFrame(stability_rows)
    topk = pd.DataFrame(topk_rows)
    consensus = (stability.assign(AbsMedianCoefficient=stability.MedianCoefficient.abs())
                 .groupby(['Threshold', 'Descriptor'], as_index=False)
                 .agg(MeanSelectionFrequency=('SelectionFrequency', 'mean'),
                      MinSelectionFrequency=('SelectionFrequency', 'min'),
                      MeanAbsMedianCoefficient=('AbsMedianCoefficient', 'mean'))
                 .sort_values(['Threshold', 'MeanSelectionFrequency', 'MeanAbsMedianCoefficient'],
                              ascending=[True, False, False]))
    consensus['Protocol'] = PROTOCOL
    for frame, name in zip((stability, topk, consensus), TABLES):
        frame.to_csv(output_dir / name, index=False)
    pd.DataFrame(audit.rows).to_csv(audit_dir / 'training_fit_audit.csv', index=False)
    pd.DataFrame(fold_rows).to_csv(audit_dir / 'nested_fold_selectors.csv', index=False)
    pd.DataFrame(tuning_rows).to_csv(audit_dir / 'C_tuning_scores.csv', index=False)
    pd.DataFrame(predictions).to_csv(audit_dir / 'test_predictions.csv', index=False)
    unchanged = all(file_hash(ROOT / path) == value for path, value in original_hashes.items())
    if not unchanged:
        raise AssertionError('Original inputs changed during analysis.')
    metadata = {'protocol': PROTOCOL, 'seed': SEED, 'alpha_grid': ALPHA_GRID, 'l1_grid': L1_GRID,
                'C_grid': C_GRID, 'topk_grid': TOPK_GRID, 'stability_repeats': STABILITY_REPEATS,
                'n_inner_validation_folds': N_FOLDS, 'n_selector_tuning_folds': N_FOLDS,
                'selectors_fitted': len(fold_rows) + len(study.splits) * len(THRESHOLDS),
                'all_topk_values_reported_without_test_selection': True,
                'inner_cv_scores_are_selection_scores_not_unbiased_performance': True,
                'fit_audit_events': len(audit.rows), 'forbidden_training_row_events': 0,
                'convergence_warnings': audit.convergence_warnings,
                'original_input_sha256': original_hashes, 'original_inputs_unchanged': unchanged,
                'script_sha256': file_hash(Path(__file__)),
                'versions': {'python': platform.python_version(), 'numpy': np.__version__,
                             'pandas': pd.__version__, 'scipy': scipy.__version__, 'sklearn': sklearn.__version__},
                'output_sha256': {name: file_hash(output_dir / name) for name in TABLES},
                'elapsed_seconds': time.perf_counter() - started}
    (audit_dir / 'protocol_and_verification.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    print(topk[['Split', 'Threshold', 'TopK', 'TestAUROC', 'BestC']].round(6).to_string(index=False))
    return metadata


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUT)
    parser.add_argument('--audit-dir', type=Path)
    args = parser.parse_args()
    run_analysis(args.output_dir, args.audit_dir)
