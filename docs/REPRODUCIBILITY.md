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
   also checks that Supplementary Tables S1-S20 are complete and in citation
   order, and reads the figure alignment and collision reports.

`docs/RESULT_FILE_MAP.md` connects the reported results to their files.
`results/final/EXPORT_PROVENANCE.json` records inputs and transformations for
the regenerated tables and figure data. The descriptor fitting ledger is
`output/descriptor_audit/training_fit_audit.csv`.
Restored model checksums, training structures and restoration comparisons are
in `output/unseen_identity_case_study/model_restore/`.

## End-to-end rerun

`python code/run_all_analysis.py --clean` deletes `results/` and reruns
scripts `00`-`11` from the raw data and split files, including the molecular
kernels. Run it in a separate copy of the repository. A full rerun in the
pinned environment and its comparison with the supplied `results/` files are
recorded in `docs/END_TO_END_RERUN.md`.

## Case-study batches

Two inference batches (`initial8` and `added2`) belong to the same case study;
together they contain ten structures and 90 scores. Six identities have the
eleven included route-specific endpoints. Both batches and their manifests are
required to reproduce this analysis.
