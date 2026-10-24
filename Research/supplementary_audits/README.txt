Supplementary source and sensitivity analyses

Run commands from the research directory with the recorded Python environment.
The analyses use the packaged benchmark and source records.

Reconstruct the label-rule, identifier and source-evidence ledgers:
  python supplementary_audits/statistics/verify_audit_statistics.py
  python supplementary_audits/identity_sources/audit_identity_sources.py compare
  python supplementary_audits/identity_sources/build_source_evidence.py

Reproduce the model sensitivities and supplementary tables:
  python supplementary_audits/statistics/run_uniform_d_exclude.py
  python supplementary_audits/statistics/run_current_identifier_flagged_exclusion.py
  python supplementary_audits/statistics/make_statistics_si.py
  python scripts/make_crossover_figure.py
  python scripts/assemble_audit_supplement.py

PubChem responses and their request metadata are preserved under
identity_sources/responses/. The fetch_* modes explicitly access the network;
later responses are not assumed to reproduce the benchmark's collection-time
records. The upstream/ directory contains source-code evidence snapshots.
The Rothamsted report is stored under identity_sources/source_documents/;
its recorded checksum is checked before use. An optional download is explicit.

The 33-identity screening set was fixed before its sensitivity fits. A flagged
identity is not automatically a mislabeled compound. All exclusions retain
paired training and evaluation membership within each comparison. Conditional
uncertainty, source limitations and the full result grids are retained.

Tables S25a-d report exclusion-log scope, aggregation-rule comparisons,
uniform-D exclusion and the Benomyl qualifier sensitivity. Tables S26a-c
report identifier comparisons, repeated-CID groups and screening-set exclusion.
