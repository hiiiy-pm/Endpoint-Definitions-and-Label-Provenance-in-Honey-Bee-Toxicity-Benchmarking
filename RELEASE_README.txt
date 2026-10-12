IEEE Access submission research materials
Endpoint Definitions and Label Provenance in Honey-Bee Toxicity Benchmarking
Prepared for public release: 12 October 2026

Authors: Yongye Huang, Runqing Liao and Zhining Wang
School of Mathematics and Statistics, Hanshan Normal University, China
Corresponding author and materials contact: Zhining Wang
Email: mathwangzhining@163.com

PURPOSE AND STATUS
These materials support the submitted manuscript's endpoint comparisons,
label-provenance audits, training-by-evaluation crossover and sensitivity
analyses. This is a submission-materials release; it does not assert acceptance,
publication or an IEEE reproducibility badge.

The 806 research files match the four final Research ZIP attachments exactly.
The repository places the contents of Research/ at its root. This additional
release README is distribution metadata outside that frozen research snapshot.
FILE_MANIFEST.json and SHA256SUMS.txt identify the archived research files.
Historical run-specific hashes and timestamps describe those original runs;
they are not rewritten as claims about this release.

This release extends historical base commit
27f63e1c3cfc63e95c60f6c3cce3edaec05af623. The base remains available in Git history.
The additions include current identifier and source audits, uniform-aggregation
and identity-exclusion sensitivities, Figure 7 and SI S25-S26 materials.
Current manuscript Tables 1-4 and the historical Table 5 / current SI S18a-b
mapping are documented in README.md, docs/REPRODUCIBILITY.md and
supplementary_audits/FINAL_MANUSCRIPT_MAP.txt.

FILES AND SIZES
Research_Code_and_Data.zip       32,774,716 bytes
Research_Model_Outputs.zip       18,690,813 bytes
Research_Source_Audits.zip        31,704,644 bytes
Research_Figures_and_Tables.zip   12,687,927 bytes
Total compressed research data:  95,858,100 bytes (about 95.9 MB)
Total uncompressed file content: 140,591,901 bytes (about 140.6 MB)
Allow additional space for the Python environment and any regenerated outputs.

All four ZIPs have non-overlapping members under Research/. Extract all four
into the same destination, then run commands from the resulting Research/
directory. A repository checkout already places those files at the root.
RESEARCH_ARCHIVE_SHA256SUMS.txt provides the four download checksums.

The four components contain, respectively:
1. Code, frozen inputs, documentation, supplementary audits and verification.
2. Molecular representations, predictions, metrics and saved model outputs.
3. External-source reconstruction, external evaluation and unseen-case evidence.
4. Figures, figure source data, tables and figure-quality records.

PLATFORM, DEPENDENCIES AND INSTALLATION
The documented end-to-end run used Windows and Python 3.12.13; see
docs/END_TO_END_RERUN.md for its actual date, environment and scope.
The current packaging/integrity verification was also performed on Windows.
Compatibility with other operating systems has not been established by this
release operation. Exact dependency versions are listed in SOFTWARE_VERSIONS.json,
requirements-analysis.txt and requirements-figures.txt; environment.yml provides
the corresponding conda specification. Historical environment records are retained.

Example installation from the research root:
  python -m venv .venv
  .venv\Scripts\python.exe -m pip install -r requirements-analysis.txt -r requirements-figures.txt

USAGE AND EXPECTED OUTPUTS
Check the shipped files without fitting models:
  .venv\Scripts\python.exe verification/verify_files.py
Expected output: status PASS, files_checked 804, empty missing and changed lists.
The manifest and checksum-list files themselves are not among those 804 entries.

For saved-result metric checks, run in a separate working copy because report
and table outputs are regenerated:
  python verification/verify_metrics.py
  python verification/verify_rank_statistics.py
These recompute metrics, comparisons and checks from saved results. Read the
scientific limitations together with the PASS/check reports; a numerical check
does not establish population-level validation.

Figure 7 and the added SI table generation entry points are:
  python scripts/make_crossover_figure.py
  python scripts/assemble_audit_supplement.py
Their expected figure/table outputs and provenance are mapped in
supplementary_audits/FINAL_MANUSCRIPT_MAP.txt. Added-analysis instructions are
in supplementary_audits/README.txt and docs/REPRODUCIBILITY.md.

For a full analysis rerun, follow docs/REPRODUCIBILITY.md in a separate copy.
Its --clean entry point deletes that copy's results/ before rebuilding it.
No new model-training or experimental run was performed to publish this release.

LICENSE AND SOURCE ATTRIBUTION
Preserve LICENSE and NOTICE.md. The MIT license applies to the authors' original
analysis code within its stated scope. Third-party data, cached upstream code,
regulatory records and source documents retain their original terms and attribution;
this release does not relicense them as MIT. Source ledgers and provenance documents
identify their origins and the limits of the reconstruction.

JOURNAL SUPPLEMENTARY-MATERIALS ROUTE
The four Research ZIPs are also intended for separate submission as journal
supplementary materials. GitHub publication complements that route and does not
replace uploading the attachments to the IEEE Access submission system.
The release documentation follows the scope of the official guidance:
https://ieeeaccess.ieee.org/authors/submission-guidelines/
https://journals.ieeeauthorcenter.ieee.org/create-your-ieee-journal-article/prepare-supplementary-materials/
