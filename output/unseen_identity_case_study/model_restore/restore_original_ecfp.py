"""Restore exactly 9 original ECFP SVCs from original development data only.

This DOES perform 9 original-data fits. It is model restoration, not new tuning;
no external structure or endpoint file is opened. The earlier probe also fitted
the original nine models, so no claim of zero training is appropriate.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
import sys

import joblib
import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.svm import SVC
from rdkit import rdBase

sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
sys.path.insert(0, str(ROOT / 'code'))
from _common import load_study_data  # Actual original split file paths and order.
from predict_external import (cross_kernel, fingerprint_generator, fingerprints_from_smiles,
                              predict_structure_table, sha256)

CHECKS = {}


def check(name, condition):
    CHECKS[name] = bool(condition)
    if not condition:
        raise AssertionError(name)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    input_paths = [ROOT / 'code/_common.py', ROOT / 'code/01_primary_endpoint_models.py',
                   ROOT / 'data/raw/apistox.csv', ROOT / 'results/cache/kernel_ecfp.npy',
                   ROOT / 'results/cache/processed_data.csv', ROOT / 'results/cache/CACHE_METADATA.json',
                   ROOT / 'results/primary/02_test_predictions.csv',
                   ROOT / 'output/external_data_audit/apistox_structures.csv',
                   ROOT / 'requirements-analysis.txt']
    input_paths += [ROOT / f'data/official_splits/{s}_{p}.csv'
                    for s in ('random', 'maxmin', 'time') for p in ('train', 'test')]
    original_hashes = {str(p.relative_to(ROOT)): sha256(p) for p in input_paths}
    study = load_study_data()
    processed = pd.read_csv(ROOT / 'results/cache/processed_data.csv')
    reference = pd.read_csv(ROOT / 'results/primary/02_test_predictions.csv')
    audit = pd.read_csv(ROOT / 'output/external_data_audit/apistox_structures.csv').fillna('')
    cached_kernel = np.load(ROOT / 'results/cache/kernel_ecfp.npy')
    check('original_size_1035', len(study.df) == 1035)
    check('original_data_processed_same_order', study.df.SMILES.tolist() == processed.SMILES.tolist())
    check('processed_tier_identical', np.array_equal(study.tier, processed.Tier.to_numpy()))
    for threshold in (100, 11, 1):
        check(f'processed_label_identical_{threshold}', np.array_equal(study.y[threshold], processed[f'Y{threshold}']))
    all_fps, identities = fingerprints_from_smiles(study.df.SMILES)
    check('all_original_identity_resolved', identities.inchikey.ne('').all() and len(identities) == 1035)
    audit_ordered = audit.set_index('SMILES').loc[study.df.SMILES]
    for field in ('inchikey', 'connectivity', 'parent_connectivity'):
        check(f'original_identity_matches_previous_audit_{field}', identities[field].tolist() == audit_ordered[field].tolist())
    regenerated = cross_kernel(all_fps, all_fps)
    check('independent_smiles_kernel_exactly_matches_cache', np.array_equal(regenerated, cached_kernel))
    check('kernel_float32_matches_original_path', regenerated.dtype == np.float32 and cached_kernel.dtype == np.float32)
    reference_identity = pd.concat([study.df[['CID', 'name', 'CAS', 'SMILES']].reset_index(names='Index'), identities], axis=1)
    reference_identity.to_csv(OUT / 'original_structure_reference.csv', index=False, encoding='utf-8-sig')
    models, restored_rows, metrics, artifacts = [], [], [], ['original_structure_reference.csv', 'predict_external.py', 'restore_original_ecfp.py']
    for split, (train, test) in study.splits.items():
        check(f'{split}_size_and_disjoint', len(train) == 828 and len(test) == 207 and not set(train) & set(test))
        training = reference_identity.iloc[train].copy().reset_index(drop=True).reset_index(names='TrainOrder')
        training['Tier'] = study.tier[train]
        for threshold in (100, 11, 1):
            training[f'Y{threshold}'] = study.y[threshold][train]
        train_filename = f'ordered_training_{split.lower()}.csv'
        training.to_csv(OUT / train_filename, index=False, encoding='utf-8-sig')
        artifacts.append(train_filename)
        for threshold in (100, 11, 1):
            model = SVC(C=1, kernel='precomputed', class_weight='balanced')
            model.fit(cached_kernel[np.ix_(train, train)], study.y[threshold][train])
            check(f'{split}_{threshold}_classes_0_1', np.array_equal(model.classes_, [0, 1]))
            check(f'{split}_{threshold}_training_feature_width', model.n_features_in_ == 828)
            original = reference[reference.Split.eq(split) & reference.Representation.eq('ECFP') & reference.Threshold.eq(threshold)].copy()
            check(f'{split}_{threshold}_original_test_sequence', original.sort_values('TestOrder').Index.tolist() == test.tolist())
            original = original.set_index('Index').loc[test]
            check(f'{split}_{threshold}_original_test_CID', original.CID.tolist() == study.df.iloc[test].CID.tolist())
            check(f'{split}_{threshold}_original_test_Y', np.array_equal(original.Y, study.y[threshold][test]))
            cache_score = model.decision_function(cached_kernel[np.ix_(test, train)])
            independent_score = model.decision_function(regenerated[np.ix_(test, train)])
            max_error = float(np.max(np.abs(independent_score - original.Score.to_numpy())))
            check(f'{split}_{threshold}_saved_score_reproduced_1e_12', max_error < 1e-12)
            check(f'{split}_{threshold}_independent_features_same_score', np.array_equal(cache_score, independent_score))
            check(f'{split}_{threshold}_fixed_class_identical', np.array_equal(independent_score > 0, original.Score.to_numpy() > 0))
            filename = f'ecfp_{split.lower()}_{threshold}.joblib'
            joblib.dump(model, OUT / filename, compress=3)
            artifacts.append(filename)
            models.append({'split': split, 'threshold_ug_bee': threshold, 'file': filename,
                           'training_file': train_filename, 'n_train': 828,
                           'ordered_training_indices': train.tolist(),
                           'n_training_positive': int(study.y[threshold][train].sum()),
                           'n_support': model.n_support_.tolist(), 'classes': model.classes_.tolist(),
                           'model_parameters': model.get_params(),
                           'positive_label_definition': f'Original Tier >= { {100: 1, 11: 2, 1: 3}[threshold] }, nominal LD50 <= {threshold} ug/bee',
                           'max_old_test_score_abs_error': max_error})
            metrics.append({'split': split, 'threshold_ug_bee': threshold, 'n_test': 207,
                            'max_abs_score_error': max_error, 'classification_disagreements': 0,
                            'independent_feature_kernel_cache_exact': True})
            for i, index in enumerate(test):
                restored_rows.append({'split': split, 'threshold_ug_bee': threshold, 'TestOrder': i,
                                      'Index': int(index), 'CID': int(study.df.iloc[index].CID),
                                      'inchikey': identities.iloc[index].inchikey,
                                      'original_saved_score': float(original.iloc[i].Score),
                                      'restored_cache_score': float(cache_score[i]),
                                      'restored_independent_features_score': float(independent_score[i]),
                                      'absolute_error': float(abs(independent_score[i] - original.iloc[i].Score)),
                                      'original_Y': int(original.iloc[i].Y), 'restored_prediction': int(independent_score[i] > 0)})
    manifest = {
        'purpose': 'Exact restoration of original primary ECFP models, enabling structure-only inference',
        'model_restoration_fits_in_this_script': 9,
        'preliminary_feasibility_probe_original_data_fits': 9,
        'not_new_hyperparameter_training_or_tuning': True,
        'external_structures_used_in_fit': 0, 'external_labels_opened_or_used_in_fit': 0,
        'original_input_sha256': original_hashes,
        'artifact_sha256': {filename: sha256(OUT / filename) for filename in artifacts},
        'fingerprint': {'type': 'Morgan bit fingerprint', 'radius': 2, 'fpSize': 1024,
                        'includeChirality': False, 'other_options': 'RDKit GetMorganGenerator defaults',
                        'generator_info': fingerprint_generator().GetInfoString(),
                        'standardization': 'None for model features; FragmentParent/Uncharger only for exclusion identities'},
        'kernel': {'type': 'Tanimoto', 'dtype_before_SVC': 'float32', 'n_original_molecules': 1035,
                   'independent_regeneration_equals_original_cache': True,
                   'raw_regenerated_kernel_sha256': hashlib.sha256(regenerated.tobytes()).hexdigest()},
        'frozen_decision_rule': 'score > 0 => positive; positive labels reflect original nested Tier definitions, not calibrated probabilities',
        'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
                     'scipy': scipy.__version__, 'scikit_learn': sklearn.__version__, 'rdkit': rdBase.rdkitVersion,
                     'joblib': joblib.__version__},
        'models': models,
        'input_interface': {'required': ['compound_id', 'SMILES'],
                            'optional': ['name', 'inchikey', 'CAS', 'source'],
                            'all_other_columns': 'Rejected; endpoint truth must be joined after inference.'},
        'limitations': 'Restoration does not improve models or certify external endpoint/identity quality. Novelty checks use the entire original development dataset. Applicability similarity is descriptive; no threshold is selected from external outcomes.'
    }
    (OUT / 'restoration_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    # End-to-end public interface: reread model artifacts and generate fingerprints
    # from old test SMILES, without reading the kernel cache inside inference.
    interface_rows = []
    for split, (_, test) in study.splits.items():
        input_table = reference_identity.iloc[test][['SMILES', 'name', 'inchikey', 'CAS']].copy()
        input_table.insert(0, 'compound_id', [f'original_index_{i}' for i in test])
        public = predict_structure_table(input_table, OUT, selected_splits=[split])
        for threshold in (100, 11, 1):
            got = public[public.threshold_ug_bee.eq(threshold)].set_index('compound_id').loc[input_table.compound_id]
            original = reference[reference.Split.eq(split) & reference.Representation.eq('ECFP') & reference.Threshold.eq(threshold)].set_index('Index').loc[test]
            error = float(np.max(np.abs(got.score.to_numpy() - original.Score.to_numpy())))
            check(f'{split}_{threshold}_public_inference_saved_score_reproduced', error < 1e-12)
            check(f'{split}_{threshold}_public_inference_frozen_decision_rule', np.array_equal(got.predicted_label, original.Score.to_numpy() > 0))
            check(f'{split}_{threshold}_public_inference_nearest_similarity', np.array_equal(got.nearest_training_tanimoto.to_numpy(), regenerated[np.ix_(test, study.splits[split][0])].max(axis=1)))
            interface_rows.append({'split': split, 'threshold_ug_bee': threshold, 'n_test': 207,
                                   'public_interface_max_abs_error': error})
    # Guard tests confirm external outcome columns cannot silently enter inference.
    valid_input = pd.DataFrame({'compound_id': ['guard'], 'SMILES': [study.df.SMILES.iloc[0]]})
    for name, bad in [('external_label', valid_input.assign(y_11=0)),
                      ('duplicate_id', pd.concat([valid_input, valid_input], ignore_index=True)),
                      ('invalid_smiles', pd.DataFrame({'compound_id': ['x'], 'SMILES': ['INVALID']})),
                      ('mismatched_key', valid_input.assign(inchikey='AAAAAAAAAAAAAA-BBBBBBBBBB-C'))]:
        try:
            predict_structure_table(bad, OUT)
        except ValueError:
            check(f'guard_rejects_{name}', True)
        else:
            check(f'guard_rejects_{name}', False)
    try:
        predict_structure_table(valid_input, OUT, require_new_structures=True)
    except ValueError:
        check('guard_rejects_original_identity_under_require_new', True)
    else:
        check('guard_rejects_original_identity_under_require_new', False)
    check('original_input_hashes_unchanged', all(sha256(ROOT / path) == value for path, value in original_hashes.items()))
    for filename, expected in manifest['artifact_sha256'].items():
        check(f'artifact_hash_verified_{filename}', sha256(OUT / filename) == expected)
    pd.DataFrame(restored_rows).to_csv(OUT / 'restored_vs_original_test_scores.csv', index=False, encoding='utf-8-sig', float_format='%.17g')
    pd.DataFrame(metrics).to_csv(OUT / 'restoration_metric_summary.csv', index=False, encoding='utf-8-sig', float_format='%.17g')
    pd.DataFrame(interface_rows).to_csv(OUT / 'public_interface_validation.csv', index=False, encoding='utf-8-sig', float_format='%.17g')
    verification = {'round': 1, 'checks': CHECKS, 'n_checks_passed': sum(CHECKS.values()),
                    'model_restoration_fits': 9, 'preliminary_original_data_probe_fits': 9,
                    'external_data_in_training': False, 'external_labels_read': False,
                    'saved_models': len(models), 'compared_old_test_scores': len(restored_rows),
                    'independent_new_structure_feature_path_validated': True,
                    'max_score_abs_error': float(max(row['max_abs_score_error'] for row in metrics)),
                    'versions': manifest['versions'], 'input_sha256': original_hashes,
                    'manifest_sha256': sha256(OUT / 'restoration_manifest.json')}
    (OUT / 'round1_restoration_verification.json').write_text(json.dumps(verification, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k, v in verification.items() if k not in ('checks', 'input_sha256')}, indent=2))


if __name__ == '__main__':
    main()
