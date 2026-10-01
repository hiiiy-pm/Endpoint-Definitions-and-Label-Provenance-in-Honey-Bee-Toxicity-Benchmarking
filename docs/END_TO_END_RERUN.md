# End-to-end rerun

Environment: Python 3.12.13 with the pinned packages of
`requirements-analysis.txt` and `requirements-figures.txt`
(numpy 2.5.3, pandas 3.0.1, scipy 1.18.1,
scikit-learn 1.8.0, RDKit 2025.09.4, statsmodels 0.14.6) on Windows 11.

Command, run in a separate copy of the repository:

```sh
python code/run_all_analysis.py --clean
```

`--clean` deletes `results/`, so the molecular kernels, descriptors and
scaffold keys are rebuilt from the raw SMILES before scripts `00`-`11` run.
All twelve scripts completed without error (about 3.5 min of wall time).

## Comparison with the supplied `results/`

| Outcome | Files |
|---|---|
| Byte-identical | 39 of 81 |
| Different bytes, identical values | 31 (line endings, float text formatting, local paths in `cache/CACHE_METADATA.json`, and input hashes in `final/EXPORT_PROVENANCE.json`) |
| Floating-point differences only | 11, all at or below 4.2e-13 |

For every differing CSV the shape, column names, text columns and missing-value
pattern were identical. The largest numeric difference, 4.2e-13, is in SVM
decision scores; all AUROC, contrast, interval, permutation and learning-curve
values reported in the paper are unchanged. The tier-boundary outputs of
scripts `10` and `11` were byte-identical.
