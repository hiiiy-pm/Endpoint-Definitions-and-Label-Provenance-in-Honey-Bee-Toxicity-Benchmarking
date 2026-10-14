# Reproducing the analyses

Run the commands below from the research directory with Python 3.12 and the
dependencies in `requirements-analysis.txt` and `requirements-figures.txt`.
The commands rebuild scientific outputs, so use a separate copy when retaining
the supplied results for comparison. `docs/RESULT_FILE_MAP.md` maps outputs to
the article, and `docs/DATA_PROVENANCE.md` identifies the source inputs.

## Benchmark analyses

```sh
python code/run_all_analysis.py --clean
```

`--clean` removes `results/` in the working copy and runs analyses `00`–`14`
from the raw ApisTox data and six official split files. The runner rebuilds
molecular kernels, fits the specified models, regenerates statistical outputs,
validates caches against their inputs and describes source-cohort composition.
It uses the packaged external-source records; it does not refresh databases
or re-adjudicate toxicological studies.

For the label-provenance analysis alone, run the following in order using the
supplied primary scores and bootstrap draws:

```sh
python scripts/prepare_cache.py
python code/10_tier_boundary_diagnostics.py
python code/11_qualifier_propagation_audit.py
python code/12_label_rule_audit.py
python code/13_corrected_label_crossover.py
python code/14_record_case_audit.py
```

The first command recomputes molecular representations from raw SMILES and
validates the supplied cache in the active environment without fitting models
or replacing cache arrays. Run it before reusing the cache after a Python
patch-version or dependency-version change. Environment matching remains exact;
do not manually edit cache JSON files to bypass that check.

Script `10` fits the 15 Tier1-excluded models and reuses saved primary scores
for paired comparisons. Script `11` reconstructs the 441 ECOTOX-derived
benchmark annotations from the cached upstream records. Script `12` separates
qualifier handling from aggregation. Script `13` uses identical retained
training and test identities within each crossover. Script `14` preserves
row-level source provenance and separates processing defects from unresolved
cross-study disagreements. `keep` and `exclude` retain nested endpoints;
`negative` is a non-nested extreme stress test.

## Regulatory-source reconstruction

The following offline workflow uses the packaged source workbooks, query
responses and existing primary predictions:

```sh
python -X utf8 scripts/rebuild_external_sources.py
```

It performs these five scientific steps with network requests blocked:

```sh
python -X utf8 output/external_data_audit/extract_tables.py
python -X utf8 output/external_data_audit/align_structures.py --offline --rebuild-structures
python -X utf8 output/external_data_audit/audit_external_data.py
python -X utf8 output/external_evaluation/build_frozen_members.py
python -X utf8 output/external_evaluation/frozen_evaluation.py
```

The membership builder derives the eligible cohort from source records and
official split identities without using predictions. The evaluator checks
that reconstruction before calculating metrics and paired intervals. An absent
or retryable PubChem response causes the offline run to fail instead of
changing the source evidence. UTF-8 mode handles dose-unit characters on
Windows. No models are fitted in this workflow.

Optional `--compare-to /path/to/reference/Research` compares rebuilt CSV/JSON/NPZ
outputs with a separate supplied-result copy. Run logs and comparison reports
are generated under `results/reproducibility/external_rebuild/`; these are
execution outputs, not required inputs. See `output/external_data_audit/README.md`
for the cache and source-matching rules.

## Identifier and aggregation sensitivities

```sh
python supplementary_audits/statistics/verify_audit_statistics.py
python supplementary_audits/identity_sources/audit_identity_sources.py compare
python supplementary_audits/identity_sources/build_source_evidence.py
python supplementary_audits/statistics/run_uniform_d_exclude.py
python supplementary_audits/statistics/run_current_identifier_flagged_exclusion.py
python supplementary_audits/statistics/make_statistics_si.py
```

The first three commands reconstruct the label-rule, identifier and source
ledgers from the packaged inputs. The next two fit the specified sensitivities.
The 33-identity screening set is a fixed scientific input, not a set inferred
anew from future database responses. The Benomyl qualifier interpretation is
a conditional sensitivity, not a confirmed same-study correction. The methods
and source limits are in `supplementary_audits/identity_sources/README.md`.

## Table and figure generation

With the numerical outputs present, run:

```sh
python code/09_export_submission_tables.py
python scripts/make_supp_tables.py
python scripts/make_external_tables.py
python scripts/make_source_holdout_table.py
python scripts/make_tier_boundary_tables.py
python scripts/make_audit_tables.py
python scripts/assemble_supplementary_tables.py
python scripts/make_main_tables.py
python scripts/generate_manuscript_figures.py --publication-layout
python scripts/generate_manuscript_figures.py --only FigS1_full_performance FigS2_model_decomposition FigS3_structural_diagnostics
python scripts/make_label_source_figure.py
python scripts/make_crossover_figure.py
python scripts/assemble_audit_supplement.py
```

The nesting diagnostics used in Table S23b are calculated by the table
generator from the benchmark labels, scenario definitions, official splits
and crossover predictions. A separately prepared training-nesting JSON is
not an input. The descriptor fitting ledger records the training-fold
boundaries used in feature selection and model fitting.

## Unseen-identity case study

The two inference batches (`initial8` and `added2`) together contain ten
structures and 90 model scores. Six identities have eleven eligible
route-specific endpoints. Both batches and their scientific input manifests
are needed for the analysis. The two batches are parts of one case study.

```sh
python output/unseen_identity_case_study/model_restore/predict_external.py --input output/unseen_identity_case_study/predictions/initial8/structure_only_input.csv --output predictions.csv --require-new-structures
python output/unseen_identity_case_study/evaluate_new_compound_pilot.py
```

The first command demonstrates inference without fitting. The second
recomputes the reported evaluation using both supplied prediction batches
and the endpoint ledger. To regenerate the model objects, run
`python output/unseen_identity_case_study/model_restore/restore_original_ecfp.py`;
this fits the nine specified models on their original training data.

## Manuscript compilation

The separate `LaTeX_Source.zip` contains the document sources and required
assets. From the extracted source directory, `python build.py --source-only`
compiles them. `python build.py /path/to/Research --regenerate` also regenerates
the assets from research outputs. LaTeX and BibTeX must be on PATH, or supplied
through the `PDFLATEX` and `BIBTEX` environment variables.
