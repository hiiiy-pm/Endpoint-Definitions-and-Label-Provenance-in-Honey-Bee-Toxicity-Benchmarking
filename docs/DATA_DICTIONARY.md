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
