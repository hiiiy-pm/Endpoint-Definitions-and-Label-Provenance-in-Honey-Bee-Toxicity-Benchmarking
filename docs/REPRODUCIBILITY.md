# Reproduction and verification

Run commands from the repository root using the recorded dependencies
(Python 3.12; `requirements-analysis.txt` and `requirements-figures.txt`).

1. `python verification/verify_files.py` checks every file against
   `FILE_MANIFEST.json`.
2. `python verification/verify_round1.py` recomputes saved point metrics,
   tier-pair AUROCs and the Eq. (4) identity, the Tier1-exclusion AUROCs and the
   exact decomposition, source-label disagreement counts by tier, the qualifier
   tracing totals, and regenerates the rows of the main-text Tables 3-5 and the
   compact descriptor table from the saved results.
3. `python verification/verify_round2.py` independently computes rank-based
   metrics (including the Tier1-excluded refits), bootstrap quantiles, censored
   endpoint labels, identity exclusion and case-study confusion matrices. It
   also checks that Supplementary Tables S1-S24 are complete and in numerical
   order, and reads the figure alignment and collision reports.

`docs/RESULT_FILE_MAP.md` connects the reported results to their files.
`results/final/EXPORT_PROVENANCE.json` records inputs and transformations for
the regenerated tables and figure data. The descriptor fitting ledger is
`output/descriptor_audit/training_fit_audit.csv`.
Restored model checksums, training structures and restoration comparisons are
in `output/unseen_identity_case_study/model_restore/`.

## End-to-end rerun

`python code/run_all_analysis.py --clean` deletes `results/` and reruns
scripts `00`-`14` from the raw data and split files, including the molecular
kernels. Run it in a separate copy of the repository. A full rerun in the
pinned environment and its comparison with the supplied `results/` files are
recorded in `docs/END_TO_END_RERUN.md`.

## Revision verification and source reconstruction

`python verification/verify_revision.py` checks cache lineage and negative
failure paths, rebuilds all molecular representations from raw SMILES, refits
all 45 primary models and both training-label models in the 225 retained-cohort
(`exclude`) crossover cells. It checks all exported crossover labels, AUROCs and
gap calculations. It writes `results/reproducibility/`; it is a stronger, fitting-based
check than the saved-score rounds above. `python
verification/check_crossover_nesting.py` audits nesting of exported test labels
and agreement of cohort/label identity across representations.
`python verification/check_training_label_nesting.py` separately reconstructs
all 45 training-label cohorts, with no fitting. Its results distinguish the
nested `keep`/`exclude` corrections from the declared non-nested `negative`
stress test and are saved in `results/reproducibility/`.

External-source preparation is a separate workflow using the frozen files in
`data/external_raw/`, retained query responses, and existing primary predictions.
Run it in a separate copy because it regenerates that copy's external tables.
The complete offline entry point is:

```sh
python -X utf8 verification/rebuild_external_sources.py
```

To compare the fresh outputs with an untouched release, add
`--compare-to /path/to/untouched/repo`. The runner blocks network requests,
rebuilds all workbook/structure caches, preserves frozen inputs, and writes
per-stage logs plus CSV/JSON/NPZ comparisons to
`results/reproducibility/external_rebuild/`. It does not fit or tune models.
The six equivalent explicit steps are:

```sh
python -X utf8 output/external_data_audit/extract_tables.py
python -X utf8 output/external_data_audit/align_structures.py --offline --rebuild-structures
python -X utf8 output/external_data_audit/audit_external_data.py
python -X utf8 output/external_data_audit/verify_audit.py
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
The 3 October 2026 isolated reconstruction completed all six stages offline and
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
`python check_manuscript.py /path/to/research-repository` checks source/data/PDF
agreement. LaTeX and BibTeX must already be on PATH (or supplied via
PDFLATEX/BIBTEX). The standalone LaTeX ZIP contains the identical source files
at its root, without an extra enclosing directory.
Figure geometry audit helpers, when available on PYTHONPATH, are invoked at
render time; the supplied `output/figure_qa/` reports are retained independently.

## Case-study batches

Two inference batches (`initial8` and `added2`) belong to the same case study;
together they contain ten structures and 90 scores. Six identities have the
eleven included route-specific endpoints. Both batches and their manifests are
required to reproduce this analysis.
