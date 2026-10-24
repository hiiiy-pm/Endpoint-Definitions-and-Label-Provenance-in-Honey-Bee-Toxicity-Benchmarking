"""Inference only for restored original ECFP models; no fitting or labels.

Required CSV fields: compound_id, SMILES.
Optional fields: name, inchikey, CAS, source.
Any other column is rejected, keeping endpoint truth out of this interface.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

DEFAULT_MODEL_DIR = Path(__file__).resolve().parent
ALLOWED_COLUMNS = {'compound_id', 'SMILES', 'name', 'inchikey', 'CAS', 'source'}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalized_name(value):
    value = re.sub(r'\s*\(Ref:.*', '', str(value), flags=re.I)
    return re.sub('[^a-z0-9]', '', value.lower())


def normalized_cas(value):
    return re.sub('[^a-z0-9]', '', str(value).lower())


def molecular_identity(smiles):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f'Invalid or empty molecular SMILES: {smiles!r}')
    key = Chem.MolToInchiKey(mol)
    if not re.fullmatch('[A-Z]{14}-[A-Z]{10}-[A-Z]', key):
        raise ValueError(f'No valid standard InChIKey for {smiles!r}')
    # A few original molecules have valid fingerprints/full keys but cannot
    # be standardized to a parent; keep the model features and flag identity.
    # Such unresolved parents are never accepted by --require-new-structures.
    try:
        parent = rdMolStandardize.FragmentParent(mol)
        parent = rdMolStandardize.Uncharger().uncharge(parent)
        Chem.RemoveStereochemistry(parent)
        parent_key = Chem.MolToInchiKey(parent)
    except Exception:
        parent_key = ''
    return mol, {'inchikey': key, 'connectivity': key.split('-')[0],
                 'parent_connectivity': parent_key.split('-')[0],
                 'parent_status': 'resolved' if parent_key else 'failed',
                 'canonical_smiles': Chem.MolToSmiles(mol), 'fragments': len(Chem.GetMolFrags(mol))}


def fingerprint_generator():
    # Exact old feature parameters; NO parent/salt standardization for features.
    return GetMorganGenerator(radius=2, fpSize=1024)


def fingerprints_from_smiles(smiles):
    gen = fingerprint_generator()
    fps, identities = [], []
    for value in smiles:
        mol, identity = molecular_identity(value)
        fps.append(gen.GetFingerprint(mol))
        identities.append(identity)
    return fps, pd.DataFrame(identities)


def cross_kernel(query_fps, training_fps):
    # The original kernel cache was float32 before SVC converted to float64.
    # Keeping this rounding step is necessary for exact original score recovery.
    return np.asarray([DataStructs.BulkTanimotoSimilarity(fp, training_fps)
                       for fp in query_fps], dtype=np.float32)


def verify_artifact(model_dir, manifest, filename):
    path = model_dir / filename
    expected = manifest['artifact_sha256'][filename]
    if sha256(path) != expected:
        raise ValueError(f'Restored artifact SHA256 mismatch: {filename}')
    return path


def predict_structure_table(structures, model_dir=DEFAULT_MODEL_DIR,
                            require_new_structures=False, selected_splits=None):
    """Return long-format raw scores, fixed classifications and novelty diagnostics.

    No endpoint fields, fitting, scaling, calibration or model selection occur.
    `selected_splits` is only for restoration testing; default scores all 9 models.
    """
    model_dir = Path(model_dir).resolve()
    missing = {'compound_id', 'SMILES'} - set(structures.columns)
    extra = set(structures.columns) - ALLOWED_COLUMNS
    if missing or extra:
        raise ValueError(f'Structure-only interface: missing columns={sorted(missing)}, disallowed columns={sorted(extra)}')
    structures = structures.copy().fillna('')
    if structures.empty or structures.compound_id.astype(str).str.strip().eq('').any():
        raise ValueError('At least one nonempty compound_id is required.')
    if structures.compound_id.astype(str).duplicated().any():
        raise ValueError('compound_id must be unique.')
    manifest = json.loads((model_dir / 'restoration_manifest.json').read_text(encoding='utf-8'))
    verify_artifact(model_dir, manifest, 'predict_external.py')
    reference = pd.read_csv(verify_artifact(model_dir, manifest, 'original_structure_reference.csv')).fillna('')
    query_fps, identities = fingerprints_from_smiles(structures.SMILES)
    if 'inchikey' in structures:
        for supplied, computed, cid in zip(structures.inchikey, identities.inchikey, structures.compound_id):
            if supplied and supplied != computed:
                raise ValueError(f'Supplied InChIKey disagrees with SMILES for {cid}: {supplied} != {computed}')
    identity_diagnostics = []
    for i, row in identities.iterrows():
        original = structures.iloc[i]
        matches = {key: bool(row[key]) and bool(row[key] in set(reference[key]))
                   for key in ('inchikey', 'connectivity', 'parent_connectivity')}
        name_match = (bool(normalized_name(original.get('name', ''))) and
                      normalized_name(original.get('name', '')) in set(reference.name.map(normalized_name)))
        cas_match = (bool(normalized_cas(original.get('CAS', ''))) and
                     normalized_cas(original.get('CAS', '')) in set(reference.CAS.map(normalized_cas)))
        diag = {f'matches_full_development_{key}': val for key, val in matches.items()}
        diag.update(matches_full_development_name=bool(name_match), matches_full_development_CAS=bool(cas_match))
        diag['passes_conservative_full_development_exclusion'] = bool(row.parent_connectivity) and not any(diag.values())
        identity_diagnostics.append(diag)
    if require_new_structures:
        failed = [str(structures.iloc[i].compound_id) for i, d in enumerate(identity_diagnostics)
                  if not d['passes_conservative_full_development_exclusion']]
        if failed:
            raise ValueError(f'Full development identity exclusion failed for: {failed}')
    allowed_splits = {m['split'] for m in manifest['models']}
    if selected_splits is not None and not set(selected_splits) <= allowed_splits:
        raise ValueError('Unknown split requested.')
    rows = []
    for split in ('Random', 'MaxMin', 'Time'):
        if selected_splits is not None and split not in selected_splits:
            continue
        training_name = f'ordered_training_{split.lower()}.csv'
        training = pd.read_csv(verify_artifact(model_dir, manifest, training_name)).fillna('')
        if training.TrainOrder.tolist() != list(range(828)):
            raise ValueError('Training order is incomplete or altered.')
        train_fps, train_identity = fingerprints_from_smiles(training.SMILES)
        if training.inchikey.tolist() != train_identity.inchikey.tolist():
            raise ValueError('Restored training identity mismatch.')
        kernel = cross_kernel(query_fps, train_fps)
        nearest = np.argmax(kernel, axis=1)
        for entry in [m for m in manifest['models'] if m['split'] == split]:
            # Only self-produced artifacts listed in our manifest are loaded.
            model = joblib.load(verify_artifact(model_dir, manifest, entry['file']))
            scores = model.decision_function(kernel)
            predicted = (scores > 0).astype(int)
            if not np.array_equal(predicted, model.predict(kernel)):
                raise AssertionError('Fixed decision threshold 0 disagrees with SVC.predict.')
            for i, score in enumerate(scores):
                near = training.iloc[int(nearest[i])]
                rows.append({
                    'compound_id': str(structures.iloc[i].compound_id),
                    'name': structures.iloc[i].get('name', ''), 'SMILES': structures.iloc[i].SMILES,
                    **identities.iloc[i].to_dict(), **identity_diagnostics[i],
                    'split': split, 'threshold_ug_bee': entry['threshold_ug_bee'],
                    'representation': 'ECFP', 'score': float(score),
                    'fixed_decision_threshold': 0.0, 'predicted_label': int(predicted[i]),
                    'positive_label_definition': entry['positive_label_definition'],
                    'nearest_training_tanimoto': float(kernel[i, nearest[i]]),
                    'nearest_training_Index': int(near.Index), 'nearest_training_CID': int(near.CID),
                    'nearest_training_inchikey': near.inchikey, 'nearest_training_name': near['name'],
                    'matches_training_inchikey': bool(identities.iloc[i].inchikey in set(training.inchikey)),
                    'matches_training_connectivity': bool(identities.iloc[i].connectivity in set(training.connectivity)),
                    'matches_training_parent_connectivity': bool(identities.iloc[i].parent_connectivity in set(training.parent_connectivity)),
                    'model_artifact_sha256': manifest['artifact_sha256'][entry['file']],
                    'score_is_probability': False,
                })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--model-dir', type=Path, default=DEFAULT_MODEL_DIR)
    parser.add_argument('--require-new-structures', action='store_true')
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError('Output must not overwrite the input structure table.')
    input_hash = sha256(args.input)
    data = pd.read_csv(args.input, dtype={'compound_id': str})
    output = predict_structure_table(data, args.model_dir, args.require_new_structures)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, index=False, encoding='utf-8-sig', float_format='%.17g')
    if sha256(args.input) != input_hash:
        raise AssertionError('Input structure table changed.')
    print(f'Scored {data.compound_id.nunique()} structures with 9 frozen restored models: {len(output)} rows.')


if __name__ == '__main__':
    main()
