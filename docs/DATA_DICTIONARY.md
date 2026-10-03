# Data dictionary

The frozen raw file is `data/raw/apistox.csv` and contains 1,035 unique
compound records.

| Column | Description |
|---|---|
| `name` | Compound name supplied by the source dataset. |
| `CID` | PubChem compound identifier. |
| `CAS` | CAS registry number when available. |
| `SMILES` | Canonicalized molecular structure string used for modelling. |
| `source` | Source database: PPDB, ECOTOX or BPDB. |
| `year` | Publication-year proxy used by the Time split. |
| `toxicity_type` | Exposure/record type: Contact, Oral or Other. |
| `herbicide` | Binary agrochemical-use flag. |
| `fungicide` | Binary agrochemical-use flag. |
| `insecticide` | Binary agrochemical-use flag. |
| `other_agrochemical` | Binary flag for other agrochemical uses. |
| `label` | Existing EPA binary toxicity label. |
| `ppdb_level` | Existing PPDB ordinal toxicity level used in Tier mapping. |

The six frozen split files are under `data/official_splits/`. Each split has
828 training and 207 test compounds and is joined to the raw data by `SMILES`.

## Metrics and revised label audit

Legacy result-column `AUPRC` means scikit-learn **average precision (AP)**,
the step-function precision–recall summary, not trapezoidal PR-AUC. AUROC
contrasts throughout are **severe minus broad** (11−100 or 1−100).

`results/label_audit/02_compound_labels_by_rule.csv` contains one row per
compound/filter/route definition/threshold/rule. `Index` is the zero-based
row in the immutable raw benchmark; CAS is checked against that row.
Labels are `1` positive (at or below the threshold), `0` negative (above),
`-1` unresolved and `-2` no eligible record. A/B ignore operators and use
median/consensus; C/D retain operators and use consensus/interval median.
The `fixed` route uses the historical determining group; `eligible` and
`all` evaluate the minimum over retained routes, without claiming a uniquely
identified corrected determining route.

The crossover predicts each threshold separately. `keep` retains original
labels where unresolved; `exclude` drops any compound unresolved at any of
the three thresholds from both training-label versions and all evaluations;
`negative` sets unresolved labels negative at 100 only. The last policy can
violate nesting and is reported only as an extreme, non-nested stress test.
`09_crossover_predictions.csv.gz` supplies row identities, original/corrected
labels and both trained score vectors. Membership SHA256 values identify
ordered train/test index arrays encoded as little-endian int64. No label in
these sensitivity outputs overwrites `data/raw/apistox.csv`.
