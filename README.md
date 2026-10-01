# Endpoint Definitions and Label Provenance in Honey-Bee Toxicity Benchmarking

Data, analysis code and results for the study by Yongye Huang, Runqing Liao and
Zhining Wang (corresponding author), Hanshan Normal University.

The study compares three nested toxicity endpoints (100, 11 and 1 micrograms
per bee) on the same 1,035 ApisTox compounds, official splits and learners. The
broad 100-microgram endpoint is the least discriminable endpoint (ECFP MaxMin
AUROC 0.592, 0.792 and 0.795). Tier-pair decomposition and Tier1 exclusion
locate most of this deficit at the boundary between compounds annotated at
(11, 100] and above 100 micrograms per bee. Route-matched regulatory records
contradict most (11, 100] annotations at 100 micrograms per bee, and a
replication of the upstream ApisTox ECOTOX labelling shows that inequality
qualifiers of limit-test results (for example `>100`) were dropped, so these
values entered the labels as exactly 100.

## Repository structure

| Location | Content |
|---|---|
| `data/raw/` | ApisTox benchmark data |
| `data/official_splits/` | Official Random, MaxMin and Time train/test membership |
| `data/external_raw/` | OpenFoodTox, EPA retrospective and ApisTox source files with input checksums |
| `code/` | Numerical analyses `00`-`11` and `run_all_analysis.py` |
| `results/` | Predictions, metrics, bootstrap resamples and molecular representations |
| `source_data/` | Data behind every figure |
| `figures/` | Main figures (Fig1-Fig6) and supplementary figures (FigS1-FigS4) |
| `output/external_data_audit/` | Regulatory observations, identity alignment, label concordance and cached upstream ApisTox files |
| `output/external_evaluation/` | Saved-score evaluation against regulatory labels |
| `output/unseen_identity_case_study/` | Curated unseen compounds, restored models, predictions and evaluation |
| `output/descriptor_audit/` | Fitting ledger of the nested descriptor analysis |
| `output/main_tables/` | Data rows of the main-text Tables 3-5 |
| `output/latex_tables/` | Supplementary Tables S1-S20 generated from the results |
| `output/figure_qa/` | Panel-alignment and text-collision reports for all figures |
| `scripts/` | Table and figure generation |
| `verification/` | File checksums and independent numerical checks |
| `docs/` | Data dictionary, provenance, result map and reproduction notes |

## Environment

Python 3.12 with the pinned packages:

```sh
python -m venv .venv
# Activate .venv using the command for your shell.
python -m pip install -r requirements-analysis.txt -r requirements-figures.txt
```

`environment.yml` gives the equivalent conda environment and
`SOFTWARE_VERSIONS.json` the recorded package versions and random seed.

## Check the supplied results

Run from the repository root:

```sh
python verification/verify_files.py
python verification/verify_round1.py
python verification/verify_round2.py
```

The first command checks every file against `FILE_MANIFEST.json`
(`SHA256SUMS.txt` gives the same checksums for `sha256sum -c`). The two
numerical checks recompute the saved metrics, tier-pair AUROCs, the Tier1
decomposition, concordance counts, case-study confusion matrices and table
rows from the saved predictions. The second check uses rank-pair computations
without scikit-learn metric functions. Neither refits the original models.
Reports are written to `verification/results/`.

## Reproduce the analyses

`python code/run_all_analysis.py --clean` reruns all numerical analyses
(`00`-`11`) from the raw data and official splits, including the molecular
kernels; run it in a separate copy of the repository. The tier-boundary and
qualifier analyses can also be run alone:

```sh
python code/10_tier_boundary_diagnostics.py
python code/11_qualifier_propagation_audit.py
```

`10` reuses the saved held-out scores and the primary bootstrap draws and fits
only the 15 Tier1-excluded models; `11` replicates the upstream ECOTOX labelling
from `output/external_data_audit/cache/upstream_ecotox.csv`, checks it against
all 441 ECOTOX-derived benchmark labels and then re-reads the inequality
operators. See `docs/REPRODUCIBILITY.md` and `docs/END_TO_END_RERUN.md`.

## Reproduce tables and figures

```sh
python code/09_export_submission_tables.py
python scripts/make_supp_tables.py
python scripts/make_external_tables.py
python scripts/make_source_holdout_table.py
python scripts/make_tier_boundary_tables.py
python scripts/assemble_supplementary_tables.py
python scripts/make_main_tables.py
python scripts/generate_manuscript_figures.py
python scripts/make_label_source_figure.py
```

`docs/RESULT_FILE_MAP.md` lists the file behind each reported result.

## Unseen-identity models

The restored models can be applied without fitting:

```sh
python output/unseen_identity_case_study/model_restore/predict_external.py --input output/unseen_identity_case_study/predictions/initial8/structure_only_input.csv --output predictions.csv --require-new-structures
```

## Source data

- ApisTox: https://doi.org/10.5281/zenodo.13350981
- ApisTox dataset-creation code: https://github.com/j-adamczyk/ApisTox_dataset
- OpenFoodTox 3.0: https://doi.org/10.5281/zenodo.19388272
- EPA retrospective tables: https://doi.org/10.1371/journal.pone.0265962

Source identifiers, exposure routes, units and inequality qualifiers are
retained in the curated records. The EPA data set is labelled EPA in the paper;
`PLOS_resolved_contact` is its key in the code.

## License and citation

The code written for this study is released under the MIT License (`LICENSE`).
Third-party data and source documents keep their original terms and
attribution (`NOTICE.md`). `CITATION.cff` gives the citation metadata.
