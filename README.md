# Endpoint Definitions and Label Provenance in Honey-Bee Toxicity Benchmarking

Data, analysis code and results for the study by Yongye Huang, Runqing Liao and
Zhining Wang (corresponding author), Hanshan Normal University.

The study compares three nested toxicity endpoints (100, 11 and 1 micrograms
per bee) on the same 1,035 ApisTox compounds, official splits and learners. Across the five primary representations and three official splits, the
broad 100-microgram endpoint had the lowest AUROC under the specified SVM
protocols (ECFP MaxMin AUROCs 0.592, 0.792 and 0.795). Tier-pair decomposition and Tier1 exclusion
locate most of this deficit at the boundary between compounds annotated at
(11, 100] and above 100 micrograms per bee. The selected route-matched source
cohorts with determinate labels at all three thresholds showed substantial
disagreement with the released Tier1 annotations at 100 micrograms per bee.
These selected comparisons do not establish a benchmark-wide error rate; see
main Table 4 and Supplementary Tables S15-S16 and S24c-d for the cohort scope.
A replication of the upstream ApisTox ECOTOX labelling shows that inequality
qualifiers of limit-test results (for example `>100`) were dropped, so these
values entered the labels as exactly 100. A factorial audit separates qualifier
handling from aggregation. Correcting identifiable labels in a matched-member
training-by-evaluation crossover does not show a consistent reduction of the
endpoint gap. The confirmed processing defect is not thereby established as
the cause of that gap; incomplete source records and conditional uncertainty
limit causal inference.

The endpoint analysis followed exploratory modelling on these same official
splits and is not a preregistered confirmatory study. See
`docs/ANALYSIS_HISTORY.md` for preserved historical evidence and its limits.

This archive supplies the complete analyses and evidence supporting the
accompanying article. The [final study release](https://github.com/hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking/releases/tag/ieee-access-submission-final-2026-10-12)
provides the code, data, model outputs, source records, figures and tables
matching the four Research Materials attachments.

## Supplementary audits

The historical benchmark inputs, official split files and reported primary
outputs are retained. Additional analyses are under `supplementary_audits/`;
they do not replace the released benchmark labels.

The extension includes a current PubChem identifier-consistency screen,
record-level checks of the upstream exclusion log, a qualifier-direction
counterfactual, uniform-aggregation label checks and ECFP sensitivity analyses.
Current database responses are dated snapshots and are not assumed to reproduce
the original data-collection responses. A flagged identity is a review candidate,
not automatically a mislabeled compound.

The additional scripts and outputs are included in the final study release.
The accompanying manuscript source includes the audit supplement and
crossover figure.
Instructions for the added analyses are in `supplementary_audits/README.txt`.

## Research Materials attachments

This is the complete research archive accompanying the submitted manuscript.
The four ZIP attachments extract into the same `Research/` directory, with
non-overlapping files. Extract all four into one destination.

| Attachment | Contents |
|---|---|
| `Research_Code_and_Data.zip` | Code, source inputs, documentation and supplementary scientific analyses |
| `Research_Model_Outputs.zip` | Molecular representations, predictions, metrics and model-analysis outputs |
| `Research_Source_Audits.zip` | External-source reconstruction, external evaluation and unseen-identity case materials |
| `Research_Figures_and_Tables.zip` | Current figures, figure source data and tables |

`FILE_MANIFEST.json` and `SHA256SUMS.txt` identify this exact archive by file
content. Upstream identifiers and recorded analysis history retain their
original scientific meaning. The final study release contains the complete
materials listed in this manifest.

## Repository structure

| Location | Content |
|---|---|
| `data/raw/` | ApisTox benchmark data |
| `data/official_splits/` | Official Random, MaxMin and Time train/test membership |
| `data/external_raw/` | OpenFoodTox, EPA retrospective and ApisTox source files with input checksums |
| `code/` | Numerical analyses `00`-`14` and `run_all_analysis.py` |
| `results/` | Predictions, metrics, bootstrap resamples and molecular representations |
| `source_data/` | Data behind every figure |
| `figures/` | Main figures (Fig1-Fig7) and supplementary figures (FigS1-FigS4) |
| `output/external_data_audit/` | Regulatory observations, identity alignment, label concordance and cached upstream ApisTox files |
| `output/external_evaluation/` | Saved-score evaluation against regulatory labels |
| `output/unseen_identity_case_study/` | Curated unseen compounds, restored models, predictions and evaluation |
| `output/descriptor_audit/` | Fitting ledger of the nested descriptor analysis |
| `output/main_tables/` | Retained base-release table rows; see the current main/SI mapping below |
| `output/latex_tables/` | Supplementary Tables S1-S24 and the added S25-S26 audit tables |
| `scripts/` | Table and figure generation, source reconstruction and cohort description |
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

## File integrity and saved results

`FILE_MANIFEST.json` and `SHA256SUMS.txt` record the delivered files and their
SHA-256 hashes. The saved predictions, metrics, bootstrap draws and figure
source data are identified in `docs/RESULT_FILE_MAP.md`.
Historical numerical rerun and cache-validation evidence is retained in
`results/reproducibility/`, `docs/END_TO_END_RERUN.md` and
`docs/CACHE_VALIDATION.md`.

## Reproduce the analyses

`python code/run_all_analysis.py --clean` reruns all numerical analyses
(`00`-`14`) from the raw data and official splits, including the molecular
kernels; run it in a separate copy of the repository. The tier-boundary and
qualifier analyses can also be run alone:

```sh
python code/10_tier_boundary_diagnostics.py
python code/11_qualifier_propagation_audit.py
python code/12_label_rule_audit.py
python code/13_corrected_label_crossover.py
python code/14_record_case_audit.py
```

`10` reuses the saved held-out scores and the primary bootstrap draws and fits
only the 15 Tier1-excluded models; `11` replicates the upstream ECOTOX labelling
from `output/external_data_audit/cache/upstream_ecotox.csv`, checks it against
all 441 ECOTOX-derived benchmark labels and then re-reads the inequality
operators. See `docs/REPRODUCIBILITY.md` and `docs/END_TO_END_RERUN.md`.
`12` separates qualifier and aggregation rules at all three cut-offs;
`13` uses identical retained training/test identities within each crossover.
The `negative` unresolved-label policy is an extreme stress test that can
violate threshold nesting; it is not a feasible corrected nested benchmark.
`14` preserves raw-row provenance and distinguishes record-level extraction
errors from unresolved cross-study differences. Input/model identity checks
are described in `docs/CACHE_VALIDATION.md`.

## Reproduce tables and figures

```sh
python code/09_export_submission_tables.py
python scripts/make_supp_tables.py
python scripts/make_external_tables.py
python scripts/make_source_holdout_table.py
python scripts/make_tier_boundary_tables.py
python scripts/make_audit_tables.py
python scripts/assemble_supplementary_tables.py
python scripts/make_main_tables.py
python scripts/generate_manuscript_figures.py --layout-13page
python scripts/make_label_source_figure.py
python scripts/make_crossover_figure.py
python scripts/assemble_audit_supplement.py
```

The current manuscript contains Tables 1-4: endpoint construction, primary
MaxMin contrasts, Tier1 exclusion and route-matched source agreement,
respectively. The former external-AUROC Table 5 is now represented in
Supplementary Tables S18a-b. Retained base-release exporters may generate
historical tables; see `supplementary_audits/FINAL_MANUSCRIPT_MAP.txt` for the
current main/SI mapping. `docs/RESULT_FILE_MAP.md` links the current results
to their source files and identifies historical exporter names.

## External sources and manuscript compilation

The raw-benchmark pipeline consumes the retained regulatory-source audit inputs;
it does not silently re-download or re-adjudicate every external observation.
`docs/REPRODUCIBILITY.md` describes the separate external-source reconstruction
commands and their network/cache boundaries.

The formal submission archive contains the main PDF, SI PDF and `LaTeX_Source/`.
The code/data archive is supplied separately. From `LaTeX_Source/`,
`python build.py --source-only` compiles both documents using the included
assets and an existing LaTeX installation. To regenerate assets from this
repository, use `python build.py /path/to/research-repository --regenerate`.
The standalone LaTeX ZIP has the same source files at its root.

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
