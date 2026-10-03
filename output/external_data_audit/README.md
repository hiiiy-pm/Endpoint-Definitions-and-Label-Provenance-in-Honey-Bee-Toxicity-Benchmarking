# Source identity and endpoint audit

The CSV files contain the current curated observations and source-label
comparisons used in Supplementary Table S15. In the fixed, same-route cohort,
OpenFoodTox has 149 compounds and 2, 1 and 82 disagreements at 1, 11 and 100
micrograms per bee; EPA has 103 compounds and 2, 3 and 51 disagreements.
`source_identity_quarantine.csv` retains excluded identity conflicts.

To reconstruct from `data/external_raw`, run from the repository root:

```sh
python -X utf8 output/external_data_audit/extract_tables.py
python -X utf8 output/external_data_audit/align_structures.py --offline --rebuild-structures
python -X utf8 output/external_data_audit/audit_external_data.py
python -X utf8 output/external_data_audit/verify_audit.py
python -X utf8 output/external_evaluation/build_frozen_members.py
python -X utf8 output/external_evaluation/frozen_evaluation.py
```

The first two commands rebuild disposable workbook and structure pickle
caches directly from the source workbooks. The membership builder reconstructs
the missing feasibility-stage membership dependency from source records and
official splits; the evaluator independently checks it before re-evaluating saved
predictions. Run the chain in a separate copy to preserve a released output set.
`python -X utf8 verification/rebuild_external_sources.py --compare-to /path/to/release/repo`
runs all six steps with network requests blocked and compares regenerated
scientific CSVs and every saved bootstrap draw with the release. Its logs and
machine-readable evidence are in `results/reproducibility/external_rebuild/`.

Cached PubChem query responses are consolidated in
`cache/pubchem_records.json` with their original query, URL and retrieval
metadata. Existing individual response files, if generated, take precedence.
`--offline` fails if a needed response is absent or retryable; no new lookup is
sent. Omitting that option preserves the legacy interactive lookup behaviour.
The upstream ApisTox source
files are retained under `cache/`; `fetch_provenance.py` refreshes those files
from their public sources when needed.

Full InChIKey is the primary overlap key; connectivity and standardized parent
connectivity provide additional exclusion checks. Names and CAS identifiers
support lookup and conflict review. Inequality-qualified measurements remain
intervals: for example, `>100` determines a negative label at 100, whereas
`>=100` does not. Measurements are grouped by identity and exposure route;
conflicting or indeterminate groups are flagged rather than assigned point
LD50 values. Different databases may contain the same underlying studies.
