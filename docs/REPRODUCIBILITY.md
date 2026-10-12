# Reproduction and result provenance

Run commands from the repository root using the recorded dependencies
(Python 3.12; `requirements-analysis.txt` and `requirements-figures.txt`).

`FILE_MANIFEST.json` and `SHA256SUMS.txt` identify the supplied files by
SHA-256 hash. The analysis scripts below regenerate the scientific outputs;
`docs/RESULT_FILE_MAP.md` maps those outputs to the manuscript.

The current manuscript contains Tables 1-4: endpoint construction, primary
MaxMin contrasts, Tier1 exclusion and route-matched source agreement,
respectively. The former external-AUROC Table 5 is now represented in
Supplementary Tables S18a-b. Retained base-release exporters may generate
historical tables; see `supplementary_audits/FINAL_MANUSCRIPT_MAP.txt` for the
current main/SI mapping. `docs/RESULT_FILE_MAP.md` identifies current results
and historical exporter names.
`results/final/EXPORT_PROVENANCE.json` records inputs and transformations for
the regenerated tables and figure data. The descriptor fitting ledger is
`output/descriptor_audit/training_fit_audit.csv`.
Restored model checksums, training structures and restoration comparisons are
in `output/unseen_identity_case_study/model_restore/`.

## End-to-end rerun

`python code/run_all_analysis.py --clean` deletes `results/` and reruns
scripts `00`-`14` from the raw data and split files, including the molecular
kernels. It then recomputes cache-validation metadata with
`scripts/prepare_cache.py` and the descriptive cohort-structure outputs with
`scripts/describe_cohort_selection.py`. Run it in a separate copy of the
repository. A historical full rerun in the
pinned environment and its comparison with the supplied `results/` files are
recorded in `docs/END_TO_END_RERUN.md`.

## Source reconstruction

The saved cache and numerical rerun records describe the historical validation
scope in `docs/CACHE_VALIDATION.md` and `docs/END_TO_END_RERUN.md`.
The `keep` and `exclude` corrected-label policies retain nested endpoints;
the `negative` unresolved-label policy is a non-nested extreme stress test.

External-source preparation is a separate workflow using the frozen files in
`data/external_raw/`, retained query responses, and existing primary predictions.
Run it in a separate copy because it regenerates that copy's external tables.
The complete offline entry point is:

```sh
python -X utf8 scripts/rebuild_external_sources.py
```

To compare the fresh outputs with an untouched release, add
`--compare-to /path/to/untouched/repo`. The runner blocks network requests,
rebuilds all workbook/structure caches, preserves frozen inputs, and writes
per-stage logs plus CSV/JSON/NPZ comparisons to
`results/reproducibility/external_rebuild/`. It does not fit or tune models.
The five equivalent scientific steps are:

```sh
python -X utf8 output/external_data_audit/extract_tables.py
python -X utf8 output/external_data_audit/align_structures.py --offline --rebuild-structures
python -X utf8 output/external_data_audit/audit_external_data.py
python -X utf8 output/external_evaluation/build_frozen_members.py
python -X utf8 output/external_evaluation/frozen_evaluation.py
```

The membership builder restores the previously missing
`output/frozen_external_feasibility/frozen_split_external_members.csv` dependency
from audited source records and official split identities, without reading
predictions or previous evaluation outputs. Its reconstruction manifest records
current input, generator and output hashes; it is not a recovered historical
freeze. The evaluator independently reconstructs those members again and requires
exact membership/label agreement before calculating metrics and paired intervals.
UTF-8 mode prevents Windows console encodings from failing on printed dose units.
The historical 3 October 2026 isolated reconstruction, which also included a
separate internal audit stage, completed all six stages offline and
matched 37 regenerated existing outputs, including every stored external
bootstrap draw, with maximum numerical difference zero. See
`results/reproducibility/external_rebuild/EXTERNAL_REBUILD_VERIFICATION.json`
and its per-file comparison. The working scientific tables were retained;
only the restored membership dependency and explicit protocol wording metadata
were added/refreshed after comparison.

See `output/external_data_audit/README.md` for cache and network behaviour.
An absent or retryable PubChem response fails in offline mode instead of changing
the source evidence. Mutable upstream refreshes are not part of this chain. The checked
upstream Python files match commit 57c2c05de293a0f54c119744eb3ca9cb03bc030d,
as recorded in `results/record_audit/UPSTREAM_PIN_CHECK.json`. Re-running the
offline record audit does not purport to re-open every original regulatory
study; record identifiers and missing fields remain explicit.

## Manuscript sources

The formal submission package supplies `LaTeX_Source/` alongside the main and SI
PDFs. The research code/data archive is separate. From the LaTeX source directory,
`python build.py --source-only` compiles the supplied assets; `python build.py
/path/to/research-repository --regenerate` first regenerates tables and figures.
LaTeX and BibTeX must already be on PATH (or supplied via PDFLATEX/BIBTEX).
The standalone LaTeX ZIP contains the identical source files at its root,
without an extra enclosing directory.

## Case-study batches

Two inference batches (`initial8` and `added2`) belong to the same case study;
together they contain ten structures and 90 scores. Six identities have the
eleven included route-specific endpoints. Both batches and their manifests are
required to reproduce this analysis.
