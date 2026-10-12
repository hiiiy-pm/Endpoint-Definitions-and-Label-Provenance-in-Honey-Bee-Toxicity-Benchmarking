IEEE Access final submission research materials

Endpoint Definitions and Label Provenance in Honey-Bee Toxicity Benchmarking
Yongye Huang, Runqing Liao and Zhining Wang
Final study release: https://github.com/hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking/releases/tag/ieee-access-submission-final-2026-10-12

DESCRIPTION AND PACKING LIST
The four Research archives contain the current study code, source data, saved
outputs, source records, figures and tables. They contain 710 distinct files.
They reproduce the Research materials supplied with the final submission.
Internal manuscript/code/layout checks, obsolete figure drafts and previews
are excluded. The 15 core analysis scripts 00-14, data inputs, final figures,
saved predictions, representations and restored models are unchanged.

- Research_Code_and_Data.zip: 31,608,748 bytes
- Research_Figures_and_Tables.zip: 12,656,092 bytes
- Research_Model_Outputs.zip: 18,690,813 bytes
- Research_Source_Audits.zip: 31,696,681 bytes

Research_Code_and_Data.zip contains analysis code, raw/source inputs,
documentation and supplementary scientific analyses.
Research_Model_Outputs.zip contains representations, predictions, metrics,
bootstrap draws and retained numerical rerun evidence.
Research_Source_Audits.zip contains external-source reconstruction, regulatory
reassessment, source endpoint evidence and unseen-identity case materials.
Research_Figures_and_Tables.zip contains current figures, source tables and
the tables used by the manuscript and Supplementary Information.

SETUP
1. Download all four Research ZIP files and the archive SHA-256 list.
2. Compare each ZIP with RESEARCH_ARCHIVE_SHA256SUMS.txt.
3. Extract all four archives into the same destination. Each contains
   non-overlapping paths under Research/; together they form one directory.
4. Read Research/README.md and Research/docs/REPRODUCIBILITY.md.

PLATFORM AND ENVIRONMENT
Python 3.12 and the pinned packages in requirements-analysis.txt and
requirements-figures.txt are required for numerical work. environment.yml
provides a conda specification. SOFTWARE_VERSIONS.json records package
versions and random seeds. Standard-library archive extraction works on
Windows, macOS and Linux; configure long-path support on Windows if needed.

RUN INSTRUCTIONS AND EXPECTED OUTPUTS
Work from the Research directory using the pinned environment.
The supplied predictions, metric CSVs, bootstrap arrays and figures can be
inspected directly; docs/RESULT_FILE_MAP.md maps the reported results.
python code/run_all_analysis.py --clean reruns the raw-benchmark numerical
pipeline in a separate copy, producing results/. This command fits models.
scripts/rebuild_external_sources.py reconstructs the retained external-source
analyses from bundled inputs; the documented cache/network boundary applies.
The README lists separate commands for current tables and figures, restored
model inference and supplementary scientific analyses. Reproduction is
distinguished from collection or adjudication of new source observations.

The manuscript/SI TeX sources and build script are supplied separately in the
journal submission's LaTeX_Source.zip. They are not one of these four archives.
FILE_MANIFEST.json and SHA256SUMS.txt record the exact research file contents.
Repository-only release automation is outside that manifest and the ZIPs.

LICENSE AND CONTACT
Study code: MIT License. Third-party data retain their own terms; see LICENSE
and NOTICE.md. Cite the associated article and original datasets (CITATION.cff).
Corresponding author: Zhining Wang, mathwangzhining@163.com.
No new experiments, model fits or scientific results were produced while
preparing this final submission package.
