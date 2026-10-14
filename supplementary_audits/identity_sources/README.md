# Identifier and source audit

This directory contains identifier comparisons and source-provenance evidence for the ApisTox benchmark. The 33-identity screening list is a fixed input to the reported exclusion sensitivity.

## Main outputs

| File | Content |
|---|---|
| `cid_structure_ledger_1035.csv` | All 1,035 frozen rows compared with current CID properties |
| `identity_summary.json` | Coverage, RDKit version, counts and input hash |
| `current_identifier_flagged.csv` | Frozen 33-identity list for a stated exclusion sensitivity |
| `CURRENT_IDENTIFIER_FLAGGED_FREEZE.json` | Rule, timestamp, indices and list hash, fixed before sensitivity results |
| `known_5_CID_groups_10_rows_impact_ledger.csv` | Five repeated-CID groups, CAS candidates, years, split membership and use-source dependencies |
| `flagged_33_CAS_candidate_ledger.csv` | CAS lookup results for all 33 flags; candidates are not approved historical replacements |
| `upstream_exclusion_log_745_row_audit.csv` | Actual route rejection versus overinclusive CAS-level logging |
| `benomyl_record_level_impact_ledger.csv` | Seven historical benomyl records with the proposed source interpretation separated |
| `benomyl_primary_rule_conditional_effects.csv` | Conditional label propagation without editing labels or fitting models |
| `responses/` | Raw HTTP bodies with separate URLs, timestamps, headers and SHA256 metadata |
| `upstream/` | Source code downloaded from pinned ApisTox commit `57c2c05de293a0f54c119744eb3ca9cb03bc030d` |
| `source_documents/` | Official Rothamsted PDF and hash/source metadata |

## What the checks establish

All 1,030 unique frozen CIDs were returned by 26 current PubChem property requests. The 1,035 rows comprise 875 full-key matches, 127 connectivity matches with full-key differences, four fragment-parent full-key matches, two fragment-parent connectivity matches, and 27 other structure differences. The 33 rows in the last three categories are review flags. They are not 33 confirmed historical errors. A matching CID does not establish the identity of the material tested in an assay.

The five repeated-CID groups contain three clear current mappings to different named structures: silthiofam, diflumetorim and phthalide. Their CAS candidates match the frozen SMILES. Chlormequat is a cation/salt case; metiram has no discrete current PubChem CID and must remain distinct from zineb. Metominostrobin has a geometry difference. The other four rows in those groups have full-key agreement.

The old PubChem SDQ literature endpoint returned HTTP 404 HTML for all nine tested stored/candidate CIDs. This was recorded as retrieval failure, not as evidence of no literature. Current creation dates were retained as evidence but were not used to fill or replace frozen years. Historical dates, Time membership corrections and use-flag corrections have not been established. Zero established errors in those fields means the evidence remains unresolved.

The actual upstream 11-cut-off filter rejects 536 records in 79 route groups. Its log also includes 209 records in 75 retained groups because logging selects on CAS. Both subsets are available row by row. Benomyl's proposed >10 interpretation is kept conditional: underlying trial identity is incomplete and no independent two-expert adjudication is claimed.

## Reproduction

Network retrieval uses `audit_identity_sources.py fetch_all`, `fetch_known`, and `fetch_candidates` with Python and requests. Saved successful responses are reused. `audit_identity_sources.py compare` requires RDKit; the audited version was 2025.09.4. `build_source_evidence.py` uses Python's standard library and the preserved responses. It uses the packaged Rothamsted PDF and verifies its recorded SHA256; optional retrieval is requested explicitly with `--download-missing-source`. The scripts write only inside this directory.

The 33-identity list is a separate frozen input. Re-running current network lookups in the future does not authorize changing that list or treating a later response as the original 2026-10-03 response.

## Sources and limits

- [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/rest/pug/): current CID/CAS properties, with exact request URLs in each response metadata file.
- [Metiram reference-substance record](https://pubchem.ncbi.nlm.nih.gov/compound/9006-42-2): reference SID 481109147; no discrete CID reported. Its raw record is `responses/metiram_reference_substance_481109147.body`.
- [EPA chemical list distinguishing metiram and zineb](https://www.epa.gov/sites/default/files/documents/2001ebdcnew.pdf): raw PDF body and request metadata are retained.
- [Pinned upstream code](https://github.com/j-adamczyk/ApisTox_dataset/tree/57c2c05de293a0f54c119744eb3ca9cb03bc030d): record conversion, metadata derivation and split logic.
- [Rothamsted original report](https://www.era.rothamsted.ac.uk/eradoc/OCR/ResReport1977p2-57-74): Table 2, printed page 56/PDF page 3. The archived original carries its Creative Commons Attribution 4.0 notice; citation and hash accompany the copy.

This is a computational and documentary review. It does not supply expert toxicological adjudication, infer independence from database overlap, or establish that current identifiers were identical to the 2024 provider records.
