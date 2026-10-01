# Source identity and endpoint audit

The CSV files contain the current curated observations and source-label
comparisons used in Supplementary Table S12. In the fixed, same-route cohort,
OpenFoodTox has 149 compounds and 2, 1 and 82 disagreements at 1, 11 and 100
micrograms per bee; EPA has 103 compounds and 2, 3 and 51 disagreements.
`source_identity_quarantine.csv` retains excluded identity conflicts.

To reconstruct from `data/external_raw`, run from the repository root:

```sh
python output/external_data_audit/extract_tables.py
python output/external_data_audit/align_structures.py
python output/external_data_audit/audit_external_data.py
python output/external_data_audit/verify_audit.py
```

The first two commands rebuild disposable workbook and structure pickle
caches. Cached PubChem query responses are consolidated in
`cache/pubchem_records.json` with their original query, URL and retrieval
metadata. Existing individual response files, if generated, take precedence.
New or retryable lookups may use the network. The upstream ApisTox source
files are retained under `cache/`; `fetch_provenance.py` refreshes those files
from their public sources when needed.

Full InChIKey is the primary overlap key; connectivity and standardized parent
connectivity provide additional exclusion checks. Names and CAS identifiers
support lookup and conflict review. Inequality-qualified measurements remain
intervals: for example, `>100` determines a negative label at 100, whereas
`>=100` does not. Measurements are grouped by identity and exposure route;
conflicting or indeterminate groups are flagged rather than assigned point
LD50 values. Different databases may contain the same underlying studies.
