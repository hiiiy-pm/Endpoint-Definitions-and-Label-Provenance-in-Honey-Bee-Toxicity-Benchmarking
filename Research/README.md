# Auditing label provenance and endpoint definitions in a honey-bee toxicity benchmark

Data, analysis code and results for the study by Yongye Huang, Runqing Liao,
Xiaopeng Yang and Zhining Wang, School of Mathematics and Statistics,
Hanshan Normal University.

The study compares three nested toxicity endpoints (100, 11 and 1 micrograms
per bee) on the same 1,035 ApisTox compounds, official splits and learners.
ECFP is the primary molecular representation; four other representations
provide robustness analyses. The broad endpoint has the lowest AUROC in the
15 representation-by-split settings. Source-record reconstruction identifies
a qualifier-handling defect in the upstream ECOTOX processing. Matched-cohort
training--evaluation crossovers separate the effects of revising fitting and
scoring labels. All six primary ECFP change intervals include zero, so the
tested partial revisions provide no evidence of a narrower endpoint gap.
Together, the analyses connect task composition, label provenance and the
measured consequences of source-supported revisions.

This is a retrospective analysis. The endpoint comparison followed exploratory
modelling on the same benchmark and official test splits; the primary
designations do not imply preregistration. Source disagreement is not treated
as proof of benchmark-wide label error or of independent toxicological studies.

- [Code repository](https://github.com/hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking)
- [Versioned research release](https://github.com/hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking/releases/tag/research-reproducibility-2026-10-24)

## Contents

The research release contains the complete `Research/` directory, including
all inputs, numerical analyses, saved outputs and final figure PDFs. Paths in
the table below and analysis commands are relative to `Research/`.

| Location | Contents |
|---|---|
| `data/raw/` | ApisTox benchmark records |
| `data/official_splits/` | Random, MaxMin and Time training/test membership |
| `data/external_raw/` | Regulatory-source files and their input metadata |
| `code/` | Numerical analyses and the end-to-end runner |
| `results/` | Molecular representations, labels, predictions, metrics and bootstrap distributions |
| `source_data/` | Numerical data used to draw the figures |
| `figures/` | Final numbered figure PDFs and the output directory used by figure generators |
| `output/external_data_audit/` | Source-record extraction, identity alignment and label comparisons |
| `output/external_evaluation/` | Evaluation of saved scores against regulatory-source labels |
| `output/unseen_identity_case_study/` | Curated endpoints, models and structure-only inference |
| `output/descriptor_audit/` | Training-fold fitting ledger for the descriptor analysis |
| `supplementary_audits/` | Identifier screening, source evidence and sensitivity analyses |
| `scripts/` | Source reconstruction and table/figure generation |
| `docs/` | Data dictionary, provenance, result map and reproduction instructions |

`docs/figure_table_sources.csv` maps every
current article figure and table to the underlying results and generators.

## Environment

Use Python 3.12 with the pinned dependencies:

```sh
python -m venv .venv
# Activate .venv using the command for your shell.
python -m pip install -r requirements-analysis.txt -r requirements-figures.txt
```

`environment.yml` provides a conda environment specification.
`SOFTWARE_VERSIONS.json` records the analysis packages and random seed.

Before an analysis reuses the supplied molecular-representation cache, validate
it in the active environment:

```sh
python scripts/prepare_cache.py
```

Run this again after a Python patch-version or dependency-version change. The
check recomputes representations from raw SMILES and compares them with the
supplied arrays before recording current validation evidence. It does not fit
models or replace cached arrays. Do not manually edit the cache JSON files to
make environment versions match.

## Reproduce the results

The supplied predictions, metrics, bootstrap distributions and figure data can
be inspected directly. `docs/figure_table_sources.csv` maps them to the article.
For individual analyses that reuse the supplied cache, first run the cache
validation command above. The full clean run below regenerates its own cache.
To rerun the numerical analyses from raw benchmark data and official splits:

```sh
python code/run_all_analysis.py --clean
```

Run this command in a separate copy: `--clean` removes that copy's `results/`
directory before rebuilding analyses `00` through `14`, including molecular
representations. The external-source inputs under `data/` and `output/` are
retained. Source reconstruction and the supplementary sensitivities have
separate entry points documented in `docs/REPRODUCIBILITY.md`.

The label audit distinguishes `keep`, `exclude` and `negative` unresolved-label
policies. The first two preserve nested endpoints; `negative` is an extreme
stress test that can violate nesting. Table generation computes the training
and test nesting diagnostics from the scientific inputs and outputs; it does
not require a separately prepared verification report.

## Tables and figures

The article contains main Tables 1–3, Figures 1–5, Supplementary Tables S1–S26
(including lettered subdivisions) and Figures S1–S6. The complete generation
commands and the mapping from analysis export filenames to these figure
numbers are in `docs/REPRODUCIBILITY.md`.

The eleven final figure PDFs are supplied under `figures/`. Table generators
read the machine-readable results and write editable table fragments under
`output/latex_tables/`. Their exact commands are in
`docs/REPRODUCIBILITY.md`; the scientific inputs and generators are mapped in
`docs/figure_table_sources.csv`. These operations use saved results and do
not require refitting the models.

## Structure-only inference

The nine restored ECFP models can score new structures without fitting:

```sh
python output/unseen_identity_case_study/model_restore/predict_external.py --input output/unseen_identity_case_study/predictions/initial8/structure_only_input.csv --output predictions.csv --require-new-structures
```

Input columns, identity screens and output fields are documented in
`output/unseen_identity_case_study/model_restore/README.md`.

## Source data, license and citation

- ApisTox: https://doi.org/10.5281/zenodo.13350981
- ApisTox dataset-creation code: https://github.com/j-adamczyk/ApisTox_dataset
- OpenFoodTox 3.0: https://doi.org/10.5281/zenodo.19388272
- EPA retrospective tables: https://doi.org/10.1371/journal.pone.0265962

The source ledgers retain identifiers, exposure routes, units, inequality
qualifiers and source locations. `PLOS_resolved_contact` is the implementation
key for the EPA cohort. Cached source snapshots are needed for reproducible
offline reconstruction; subsequent database responses may differ.

Code written for this study is released under the MIT License (`LICENSE`).
Third-party material retains its source terms and attribution (`NOTICE.md`).
Please cite the associated article and the original datasets; `CITATION.cff`
provides the software citation metadata.
