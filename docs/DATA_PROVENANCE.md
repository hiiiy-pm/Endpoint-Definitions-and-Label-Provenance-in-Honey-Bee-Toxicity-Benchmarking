# Data provenance and integrity

- Benchmark data set: `data/raw/apistox.csv`
- Rows: 1,035
- Columns: 13
- SHA256: `de4b160a23601a44ac5ec5b23aae95973acba060135630a6df8f73d128baea9e`
- No continuous LD50 values were reconstructed from categorical labels.
- Official Random, MaxMin and Time train/test files are included under
  `data/official_splits/`. Their compound membership matches the Zenodo release
  exactly; row order differs from the release for Random and Time, because the
  files were reconstructed rather than copied byte for byte.

## External sources

The eight supplied raw files are under data/external_raw with INPUT_MANIFEST.json. OpenFoodTox v7: https://doi.org/10.5281/zenodo.19388272; ApisTox v2: https://doi.org/10.5281/zenodo.13350981; EPA retrospective: https://doi.org/10.1371/journal.pone.0265962. The case-study endpoint ledger retains every regulatory URL, study identifier, exact source location, stage evidence and material uncertainty. No source is assumed study-independent from the benchmark solely because its database differs. Sources and their terms/attribution remain identifiable.
