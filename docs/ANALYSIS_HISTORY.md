# Analysis history and designation of primary contrasts

The endpoint comparison followed exploratory modelling on the same ApisTox
benchmark and official test splits. Earlier model-comparison reports contain
held-out AUROC/MCC results. The retained paper package includes endpoint results
alongside a later writing plan calling two contrasts primary; its README states
that the package was frozen after exploratory, mechanism-validation and final
statistical-seal rounds. These materials do not establish prospective selection
of the contrasts before outcomes were inspected.

The revised manuscript therefore uses **designated primary contrasts**, describes
the study as retrospective, and does not claim preregistration. Within-fit
training/test separation is distinct from the earlier reuse of benchmark results
for research development. In particular, a paper-ready freeze is not evidence
that all sensitivity analyses were chosen before seeing endpoint performance.

Selected pre-existing records are retained under `docs/analysis_history/` with
archive and member SHA256 hashes in `EVIDENCE_MANIFEST.json`. ZIP timestamps
indicate an apparent sequence on August 28, but are local archive metadata, not
trusted timestamps or proof of the actual first execution time. The evidence
supports disclosure of an exploratory workflow; it does not prove that no other
historical document ever existed.

The correction scenarios and record audit in scripts 12–14 were developed after
inspection of the primary endpoint results and source-label disagreements. All
reported scenarios and negative findings are retained. The current revision also
repairs the exclude-policy crossover so that both training-label versions use
the same retained training identities; its earlier exclude outputs are superseded.
