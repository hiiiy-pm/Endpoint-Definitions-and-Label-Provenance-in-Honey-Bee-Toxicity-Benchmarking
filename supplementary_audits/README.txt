Dated extension to the frozen ApisTox study

Run commands from the research directory with the recorded Python environment.
These commands use the packaged inputs and do not replace the released benchmark.

Reconstruct supplementary statistics and source evidence (no new model fitting):
  python supplementary_audits/statistics/verify_audit_statistics.py
  python supplementary_audits/identity_sources/audit_identity_sources.py compare
  python supplementary_audits/identity_sources/build_source_evidence.py --check-source-only

Reproduce the added model sensitivities:
  python supplementary_audits/statistics/run_uniform_d_exclude.py
  python supplementary_audits/statistics/run_current_identifier_flagged_exclusion.py
  python supplementary_audits/statistics/make_statistics_si.py

Regenerate the added figure and audit tables:
  python scripts/make_crossover_figure.py
  python scripts/assemble_audit_supplement.py

Current PubChem responses and their request metadata are preserved under
identity_sources/responses/. The fetch_* modes explicitly access the network;
current responses are not assumed to reconstruct the historical literature dates.
The upstream/ directory contains evidence snapshots, not independent entrypoints.
The official Rothamsted source PDF is stored under identity_sources/source_documents/
and its checksum is checked before use; optional download is explicit.

The 33-identity screening set was fixed before its sensitivity fits. A flagged
identity is not automatically a mislabeled or incorrectly represented compound.
All new exclusions retain paired training and evaluation membership within each
comparison. Conditional uncertainty and the full result grids are retained.
