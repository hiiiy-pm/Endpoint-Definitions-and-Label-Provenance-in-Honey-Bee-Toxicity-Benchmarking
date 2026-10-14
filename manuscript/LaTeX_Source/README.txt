LaTeX source package

Main source: huang_access.tex
Supplementary source: supplementary.tex
Cover-letter source: cover_letter.tex
Bibliography: references.bib

Check the TeX installation and supplied assets:
  python build.py --check-only

Build all three documents from the supplied sources:
  python build.py --source-only

The source-only build requires Python, pdfLaTeX and BibTeX. It includes the
journal class, required fonts, author photographs, and manuscript figures.
PDFLATEX and BIBTEX can specify explicit executable paths; otherwise the tools
must be on PATH. Missing TeX packages are reported by the dependency check.
Install them with the package manager of the selected TeX distribution.

The study code and data are available at:
https://github.com/hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking
Complete research materials are available from:
https://github.com/hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking/releases/latest

Extract all Research Materials archives into the same parent directory,
merging their Research/ contents. To copy generated figure and table assets:
  python build.py ../Research

To regenerate tables and figures from saved research results and compile:
  python build.py ../Research --regenerate
This requires the environment specified in Research/requirements-analysis.txt
and Research/requirements-figures.txt. See Research/README.md for analytical
reruns and source-data reconstruction.

Build outputs are huang_access.pdf, supplementary.pdf and cover_letter.pdf.
The submitted PDFs are named Manuscript.pdf, Supplementary_Information.pdf
and Cover_Letter.pdf. Compilation logs are generated in build_logs/.
