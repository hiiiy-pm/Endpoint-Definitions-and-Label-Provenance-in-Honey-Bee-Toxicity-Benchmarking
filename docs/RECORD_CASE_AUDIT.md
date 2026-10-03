# Source-record and case audit

This audit was completed on 2026-10-02 for the review requests concerning
record adjudication, study identity, source selection and experimental
metadata. It uses the supplied records without fitting models or changing
benchmark labels. Run `python code/14_record_case_audit.py` from the pinned
analysis environment. The ordinary run is offline. Add `--verify-upstream`
to repeat the three read-only upstream code comparisons.

## What the evidence establishes

The 828 retained records in the historical determining CAS/route groups of
441 ECOTOX-derived compounds can be traced to exact rows of the cached
upstream export. Replaying those same rows with and without operators
isolates a documented extraction defect. A comparison between a benchmark
label and an OpenFoodTox or EPA label instead establishes source disagreement
unless the same underlying experiment is identified. Neither chemical
identity nor route agreement establishes study identity.

The ECOTOX export contains reference numbers, authors, titles and years, but
no ECOTOX test/result identifier or MRID fields. Reference 184644, used by
many examples, is the 2020 *Pesticide Ecotoxicity Database: Bee Studies*; it
is not an individual study identifier. OpenFoodTox supplies result,
document and literature UUIDs. No MRID string was found in the supplied bee
records or literature export. EPA S1/S4 supplies MRIDs for many endpoints.
Consequently, the supplied identifiers do not support a confirmed
cross-source same-experiment match. This does not show that the sources
are independent or that shared studies do not exist.

Within all 319 aligned EPA endpoints, 304 records have 206 distinct MRIDs;
98 records repeat a previously represented MRID. The deduplication ledger
groups these report identifiers while retaining the associated chemicals,
routes and source cells. In the route-matched EPA cohort, 128 of 134
records have 128 distinct MRIDs. OpenFoodTox document/result UUIDs are
reported separately rather than counted as independent physical studies.
Four byte-content-equivalent record duplicates beyond the first occur
among the 828 ECOTOX determining rows. They are retained to reproduce the
upstream aggregation; the duplicate ledger makes them explicit.

## Deterministic illustrative cases

Strata are defined by Tier and interval-label outcome. Within each stratum,
selection prioritizes groups whose records are all adult, exactly two days
and not expressed as daily doses, then the fewest records, then ascending
numeric CAS, then benchmark index. All candidates and stratum counts are
published. These are illustrative cases, not a random sample or a new
prevalence estimate. The selection was devised retrospectively for this
audit and was not preregistered.

The table uses the historical determining route. At 1 and 100 ug/bee,
operator-aware labels use the feasible interval of the group median; at
11 they require all records to be determinate and agree. `U` means
unresolved and is not a negative label. Full names, identities, source
row hashes, reference details and labels under both aggregation rules are
in the CSV ledgers.

| Case | CAS | Source values | Stage/time/material | Labels at 1, 11, 100: original → aware |
|---|---|---|---|---|
| E01 | 76-06-2 | >100 AI ug/org | Adult; 2 d; active ingredient | 0,0,1 → 0,0,0 |
| E03 | 61-31-4 | >96.7 AI ug/org | Adult; 2 d; active ingredient | 0,0,1 → 0,0,U |
| E04 | 78-34-2 | 21.27 ug/bee | Adult; 2 d; formulation | 0,0,1 → 0,0,1 |
| E05 | 66-81-9 | >241.72 ug/bee | Adult; 2 d; formulation | 0,0,0 → 0,0,0 |
| E06 | 115970-17-7 | Three >1000 ng/org records | Stage not reported; 2, 3, 4 d; active ingredient | 1,1,1 → 0,U,U |
| E07 | 119-64-2 | >0.336 ug/bee | Adult; 2 d; formulation | 1,1,1 → U,U,U |
| E08 | 4685-14-7 | >6.04 ug/bee | Adult; 2 d; formulation | 0,1,1 → 0,U,U |

Stratum E02 sought a Tier1 group becoming determinate-negative at 100
through other right-limit values without any exactly-`>100` record. It
has zero candidates; no case was substituted. The determining-record
ledger contains 106 exactly-`>100` observations and 234 other right-limit
observations across all Tiers. These are record counts, not affected
compound counts. E06 demonstrates a severe-boundary defect in the
distributed annotations but is not a harmonized adult-48-hour example.

Four further cases are selected separately: the alphabetically first full
InChIKey among 100-ug/bee disagreements in each OFT/EPA × ECOTOX/PPDB
stratum. Every external record for each selected compound is retained.
They concern pyraflufen-ethyl, quinmerac, metrafenone and cyflufenamid.
The latter has no CAS in the aligned EPA table; its benchmark CAS is
provided in a separate field, preserving rather than concealing the
source missingness. MRIDs, OFT UUIDs, original sheet cells, operators and
matched benchmark source-row ordinals are supplied. These four cases are
source disagreements, not adjudicated same-study extraction errors.

## Source selection and metadata limits

The route-matched OFT and EPA cohorts contain 216 and 133 compounds,
respectively; 71 compounds occur in both and their union contains 278.
These overlapping sources must not be treated as independent replications.
The following benchmark-wide flow applies the source screen, full-key
match, original benchmark route match and unanimity/determinacy at all
three cutoffs. It is not the split-specific saved-score reassessment flow.

| Source screen | Full-key benchmark matches | Same route | All three labels determinate |
|---|---:|---:|---:|
| OpenFoodTox single-structure screen | 250 | 216 | 149 |
| OpenFoodTox 48-hour/material-field screen | 130 | 96 | 72 |
| EPA adult-acute records | 149 | 133 | 103 |

The full flow includes record counts and source-screen denominators.
For each final cohort, included and excluded benchmark compounds are
described by source, route, Tier and each use flag. These comparisons are
descriptive; no representativeness or selection-adjusted inference is made.

Among 263 route-matched OFT records, duration is missing in 125 and equals
48 hours in 126 (12 have other durations); material text is absent in 215;
dose basis is missing in 251, explicitly active ingredient in 11, and
another basis in one. The explicit life-stage field is absent for all
263 records. An acute-contact/oral category and an empty material field
do not independently verify adult life stage or pure active ingredient.
EPA's 134 route-matched endpoint rows identify the adult-acute category,
but S1/S4 do not supply explicit per-endpoint duration, test material or
dose basis. Those fields are therefore missing in this audit rather than
inferred from MRID availability or from an identity-quality screen.

## Provenance, arithmetic and validation

The cached `upstream_ecotox.py`, `upstream_ppdb_and_bpdb.py` and
`upstream_processing.py` each match byte-for-byte the corresponding files
at upstream commit
[`57c2c05de293a0f54c119744eb3ca9cb03bc030d`](https://github.com/j-adamczyk/ApisTox_dataset/tree/57c2c05de293a0f54c119744eb3ca9cb03bc030d).
`UPSTREAM_PIN_CHECK.json` records exact raw-file URLs and local/remote
SHA256 values. No local cache was replaced.

Every exported source record has its input-file SHA256 and original parsed
data-row ordinal (zero-based, excluding the header; not a physical text
line number). The row hash is SHA256 of UTF-8 compact JSON `[column,value]`
pairs in original column order, with missing fields represented as empty
strings. OFT workbook row numbers/UUIDs and EPA source cells/MRIDs are
retained. `INPUT_SHA256.json` covers twelve direct inputs, including the
original workbooks and existing script-12 label table.

The audit independently implements median intervals using rational
arithmetic and exact open/closed endpoints, including even-sized groups
and tied bounds. All 2,646 comparisons with script 12's fixed-route,
all-record C/D labels agree. Thus no label discrepancy from script 12's
epsilon convention was detected in this scope; this is not a proof for
arbitrary future inputs or all route/filter combinations. Approximate
`~` values are explicitly flagged and treated as points, matching the
existing analysis. Daily-dose units are flagged; their historical numeric
conversion is replayed, not endorsed as equivalent to a single dose.

A separate validation command checked all twelve input-file hashes,
828 ECOTOX row ordinals/hashes/operators, all 1,070 external row hashes and
identity/route links, 751 OFT UUID/operator/document links to the original
workbook and 319 EPA values/operators against the original S1 cells.
Independently specified expected labels for all seven cases also passed.
Results are recorded in `INDEPENDENT_VALIDATION.json`.

## Outputs

All outputs are under `results/record_audit/`:

- `01` and `02`: all ECOTOX determining rows and all 441 case candidates.
- `03`–`05`: selection rules/counts, selected cases and their complete records.
- `06`–`08`: all external rows, route-matched rows and separate cross-source cases.
- `09`–`11`: identifier inventory, EPA MRID grouping and metadata missingness.
- `12`: exact-interval transition counts at all three thresholds.
- `13`–`15`: cohort flow, included/excluded composition and exact duplicate rows.
- JSON files: summary, input hashes, pinned upstream comparison and independent validation.
