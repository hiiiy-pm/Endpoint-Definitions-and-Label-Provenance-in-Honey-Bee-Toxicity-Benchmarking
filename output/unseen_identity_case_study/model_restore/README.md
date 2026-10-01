# Original ECFP model restoration and structure-only inference

This directory restores the nine original primary ECFP SVC models: three official
splits (Random, MaxMin, Time) crossed with thresholds of 100, 11, and 1 µg/bee.
It does not introduce new hyperparameters, calibration, model selection, or
external observations into fitting. **Restoration does perform original-data
fitting:** the completed script run fitted nine models, and the preliminary
feasibility probe fitted the same nine original configurations once.

## Verified restoration

- Original split paths and training row order are read through the actual
  `code/_common.py:load_study_data`, using
  `data/official_splits/`.
- Original settings are `SVC(C=1, kernel='precomputed', class_weight='balanced')`.
- Features are RDKit Morgan bit fingerprints, radius 2, 1024 bits, all remaining
  generator settings as in the original code. Molecular inputs are not changed
  by the parent-standardization steps used for identity exclusion.
- Tanimoto values are rounded to float32 before reaching SVC, matching the old
  cache and preserving its numerical feature path.
- Regenerating all 1,035 × 1,035 similarities directly from original SMILES gives
  exactly the original kernel matrix.
- All 1,863 original ECFP test scores are reproduced with a maximum absolute
  difference of 2.220446049250313e-16. Fixed-threshold classifications agree.
- The public prediction interface independently reloads models, regenerates
  fingerprints from the old test SMILES, and reproduces those same scores.
- First-round self-tests: 135 passed. See `round1_restoration_verification.json`.
  Current independent checks are in the repository-root `verification/` directory.

## Structure-only input

Required columns: `compound_id`, `SMILES`.

Optional columns: `name`, `inchikey`, `CAS`, `source`.

All other columns are rejected. Keep experimental LD50, inequality qualifiers,
derived labels, and evaluation decisions in a separate truth table and join them
to the resulting predictions by `compound_id` only after inference. A supplied
InChIKey must agree with the submitted SMILES.

CLI usage (PowerShell):

```powershell
python -X utf8 -s `
  './output/unseen_identity_case_study/model_restore/predict_external.py' `
  --input 'STRUCTURES_ONLY.csv' --output 'PREDICTIONS.csv' --require-new-structures
```

Python interface:

```python
from predict_external import predict_structure_table
predictions = predict_structure_table(structures, require_new_structures=True)
```

`--require-new-structures` excludes identities matching **any of the full original
1,035 development compounds**, using complete InChIKey, connectivity, standardized
parent connectivity, and supplied CAS/name. Unresolved parent identity also fails
this gate. This is a conservative identity screen, not proof of independent
underlying toxicological studies; external endpoint provenance still needs review.

## Output

Each input compound receives nine rows, one for each frozen restored model.
Outputs contain the raw SVC score, the fixed `score > 0` classification, threshold,
model hash, structure identifiers, and the maximum training-set ECFP Tanimoto
similarity with the nearest compound's index, CID, name, and InChIKey. Ties select
the first molecule in the original training order. Scores are not probabilities.
No applicability-domain cutoff is selected from external outcomes.

The interface verifies hashes of the model, training identity files, complete
development reference, and inference script before use. All original input hashes
and software versions are recorded in `restoration_manifest.json`. Nine model
objects and three ordered training tables are saved locally; original data,
results, manuscript files and original code are unchanged.

## Files for review

- `restored_vs_original_test_scores.csv`: all 1,863 paired original/restored scores.
- `restoration_metric_summary.csv`: nine per-model maximum score errors.
- `public_interface_validation.csv`: independent end-to-end inference recovery.
- `round1_restoration_verification.json`: first-round checks and original hashes.
- `restoration_manifest.json`: exact model settings, original training indices,
  versions, artifact hashes, fingerprint definition and provenance.
- `ordered_training_*.csv`: ordered training structures and original labels.

Do not tune these models, choose among them using external outcomes, or refit them
with the new external compounds. Independent endpoint/identity checks remain
necessary before describing any new evaluation as external validation.
