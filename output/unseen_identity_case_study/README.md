# Unseen-identity evaluation

This analysis applies nine frozen original ECFP models to ten candidate
structures. Six chemical identities supply eleven eligible 48-hour bee
endpoints: five contact and six oral. The model objects, prediction batches,
endpoint freeze manifests and complete inclusion ledger are retained here.

Run inference without fitting from the repository root:

```sh
python output/unseen_identity_case_study/model_restore/predict_external.py --input output/unseen_identity_case_study/predictions/initial8/structure_only_input.csv --output predictions.csv --require-new-structures
```

The interface accepts `compound_id,SMILES` and optional identifier fields,
verifies model hashes and excludes structures matching the development data.
`model_restore/README.md` describes the features and output columns.

Recompute the evaluation from the frozen predictions and endpoints:

```sh
python output/unseen_identity_case_study/evaluate_new_compound_pilot.py
python output/unseen_identity_case_study/verify_evaluation_round1.py
```

The evaluator checks the endpoint and prediction manifests before joining
them. It reports route-specific metrics, fixed-decision confusion matrices,
paired contrasts and 5,000 bootstrap draws. Single-class discrimination
metrics remain undefined. Current manuscript-wide checks are in the root
`verification/` directory.

`data_curation/frozen_endpoint_ledger.csv` records the source evidence and
inclusion decisions. `evaluation/new_compound_case_predictions.csv` has 99
endpoint/model rows; these represent repeated evaluations of six chemical
identities. Both prediction batches form this single pilot. Sample size and
missing severity tiers limit the interpretation of its discrimination scores.
