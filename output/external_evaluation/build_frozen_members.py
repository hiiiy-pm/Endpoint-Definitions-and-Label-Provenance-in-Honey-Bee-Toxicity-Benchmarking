"""Rebuild the missing external feasibility membership from audited raw-source tables.

This independent grouping implementation supplies the reference table subsequently
checked by frozen_evaluation.py. It never reads frozen evaluation outputs or scores.
The new manifest records this reconstruction, not an invented historical freeze.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import json
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / 'output/external_data_audit'
OUT = ROOT / 'output/frozen_external_feasibility'
SPLITS = ('random', 'maxmin', 'time')
SOURCE_NAMES = ('OFT_single_structure_contact', 'OFT_48h_contact', 'PLOS_resolved_contact')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalized(value, field):
    value = str(value)
    if field == 'name':
        value = re.sub(r'\s*\(Ref:.*', '', value, flags=re.I)
    if field in ('CAS', 'name'):
        value = re.sub(r'[^a-z0-9]', '', value.lower())
    return value


def main():
    paths = [AUDIT / 'apistox_structures.csv', AUDIT / 'oft_acute_endpoints_aligned.csv',
             AUDIT / 'plos_adult_acute_endpoints_aligned.csv']
    paths += [ROOT / f'data/official_splits/{split}_{part}.csv'
              for split in SPLITS for part in ('train', 'test')]
    before = {path.relative_to(ROOT).as_posix(): sha(path) for path in paths}
    benchmark, oft, epa = (pd.read_csv(path).fillna('') for path in paths[:3])
    assert len(benchmark) == 1035 and benchmark.SMILES.is_unique and benchmark.inchikey.is_unique
    sources = dict(zip(SOURCE_NAMES, (
        oft[oft.single_structure_screen & oft.route.eq('Contact')],
        oft[oft['48h_no_material_flag_screen'] & oft.route.eq('Contact')],
        epa[epa.inchikey.ne('') & epa.route.eq('Contact')],
    )))
    rows, flows = [], []
    for split in SPLITS:
        train_smiles = pd.read_csv(ROOT / f'data/official_splits/{split}_train.csv').SMILES
        test_smiles = pd.read_csv(ROOT / f'data/official_splits/{split}_test.csv').SMILES
        assert train_smiles.is_unique and test_smiles.is_unique
        assert len(train_smiles) == 828 and len(test_smiles) == 207
        assert not set(train_smiles) & set(test_smiles)
        assert set(train_smiles) | set(test_smiles) == set(benchmark.SMILES)
        train = benchmark[benchmark.SMILES.isin(train_smiles)]
        test_keys = set(benchmark.loc[benchmark.SMILES.isin(test_smiles), 'inchikey'])
        for source, external in sources.items():
            matched = np.zeros(len(external), dtype=bool)
            for field in ('inchikey', 'connectivity', 'parent_connectivity', 'CAS', 'name'):
                training = {normalized(value, field) for value in train[field]} - {''}
                matched |= external[field].map(lambda value: normalized(value, field)).isin(training).to_numpy()
            blocked_keys = set(external.loc[matched, 'inchikey'])
            eligible = external[~external.inchikey.isin(blocked_keys)]
            # Any unknown/conflicting observation leaves its compound label unresolved.
            for key, group in eligible.groupby('inchikey', sort=True):
                item = {'source': source, 'split': split, 'inchikey': key}
                for threshold in (1, 11, 100):
                    values = group[f'y_{threshold}'].to_numpy(dtype=int)
                    assert set(values) <= {-1, 0, 1}
                    item[f'y_{threshold}'] = int(values.min()) if values.min() >= 0 and values.min() == values.max() else -1
                item['exact_match_to_original_test'] = key in test_keys
                rows.append(item)
            flows.append({'source': source, 'split': split, 'source_records': len(external),
                          'training_identity_matching_records': int(matched.sum()),
                          'eligible_records': len(eligible), 'eligible_keys': eligible.inchikey.nunique()})
    table = pd.DataFrame(rows)
    assert not table.duplicated(['source', 'split', 'inchikey']).any()
    assert table.groupby(['source', 'split']).ngroups == 9
    assert all(sha(ROOT / name) == value for name, value in before.items())
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / 'frozen_split_external_members.csv'
    table.to_csv(target, index=False, encoding='utf-8', lineterminator='\n')
    manifest = {'status': 'passed', 'provenance': 'Reconstructed during the 2026-10-03 revision; not a recovered historical preregistration or freeze',
                'generator': str(Path(__file__).relative_to(ROOT)).replace('\\', '/'),
                'generator_sha256': sha(Path(__file__)), 'input_sha256': before,
                'output_sha256': {target.relative_to(ROOT).as_posix(): sha(target)},
                'cohorts': len(flows), 'member_rows': len(table), 'flow': flows,
                'models_fitted': 0, 'prediction_files_read': 0,
                'identity_exclusion': 'Exclude the entire external full key if any source record matches a training full key, connectivity, parent connectivity, CAS or normalized name',
                'labels': 'All source rows within a retained key must agree on a determinate label at each cutoff'}
    (OUT / 'REBUILD_MANIFEST.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({'status': 'passed', 'cohorts': len(flows), 'member_rows': len(table)}, indent=2))


if __name__ == '__main__':
    main()
