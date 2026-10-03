"""Independently reconstruct and check all crossover training-label cohorts; no fitting."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

import verify_revision as verification

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/reproducibility'


def main():
    study = verification.c.load_study_data()
    scenarios = verification.crossover_labels(study)
    rows = []
    for scenario, correction in scenarios.items():
        eligible = np.logical_and.reduce([correction[t] >= 0 for t in verification.c.THRESHOLDS])
        for policy in ('keep', 'exclude', 'negative'):
            keep = eligible if policy == 'exclude' else np.ones(len(study.df), bool)
            labels = {}
            for threshold in verification.c.THRESHOLDS:
                replacement = 0 if policy == 'negative' and threshold == 100 else study.y[threshold]
                labels[threshold] = (correction[threshold] if policy == 'exclude'
                                     else np.where(correction[threshold] < 0, replacement, correction[threshold]))
            for split, (train, _) in study.splits.items():
                indices = train[keep[train]]
                violation_1_11 = labels[1][indices] > labels[11][indices]
                violation_11_100 = labels[11][indices] > labels[100][indices]
                violation = violation_1_11 | violation_11_100
                assert np.all(study.y[1][indices] <= study.y[11][indices])
                assert np.all(study.y[11][indices] <= study.y[100][indices])
                if policy != 'negative':
                    assert not violation.any(), (scenario, policy, split, indices[violation])
                rows.append({'Scenario': scenario, 'Unresolved': policy, 'Split': split,
                             'n_train': len(indices), 'y1_gt_y11': int(violation_1_11.sum()),
                             'y11_gt_y100': int(violation_11_100.sum()), 'violations': int(violation.sum()),
                             'indices': ';'.join(map(str, indices[violation]))})
    frame = pd.DataFrame(rows)
    assert len(frame) == 45
    inputs = [ROOT / 'data/raw/apistox.csv', ROOT / 'results/label_audit/02_compound_labels_by_rule.csv',
              ROOT / 'output/external_data_audit/apistox_structures.csv',
              ROOT / 'output/external_data_audit/same_route_label_concordance.csv',
              *sorted((ROOT / 'data/official_splits').glob('*.csv'))]
    report = {
        'status': 'passed_with_declared_non_nested_stress_tests',
        'method': 'Independent correction/cohort reconstruction from rule ledger, external concordance and official splits; no model fitting',
        'cohorts': len(frame), 'original_training_label_violations': 0,
        'keep_exclude_training_violations': int(frame[frame.Unresolved != 'negative'].violations.sum()),
        'negative_policy_training_cohorts_with_violations': int((frame[frame.Unresolved == 'negative'].violations > 0).sum()),
        'distinct_negative_policy_training_indices': sorted({int(i) for s in frame[frame.Unresolved == 'negative']['indices'] for i in s.split(';') if i}),
        'input_sha256': {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'interpretation': 'Negative-policy violations are declared stress-test semantics, not a coherent cumulative toxicity correction. Repeated identities across scenarios/splits are not independent observations.'}
    OUT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT / 'training_label_nesting_check.csv', index=False)
    (OUT / 'TRAINING_LABEL_NESTING_CHECK.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in report.items() if k not in ('input_sha256', 'distinct_negative_policy_training_indices')}, indent=2))


if __name__ == '__main__':
    main()
