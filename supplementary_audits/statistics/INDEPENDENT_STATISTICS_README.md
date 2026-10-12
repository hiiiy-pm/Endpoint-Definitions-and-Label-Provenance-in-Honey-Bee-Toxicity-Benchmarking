# Independent revision statistics

The standard-library replay passed all assertions. Frozen source data and results were not changed.

- All-D/keep and D/C/D/keep have identical three label vectors, with explicit byte hashes.
- All-D/exclude adds indices 861 and 953: training in Random/MaxMin, test in Time.
- The CAS-level exclusion log has 745 entries: 536 in deleted route groups and 209 in retained route groups.
- The Benomyl sign replay leaves all five scenarios' keep/exclude effective inputs unchanged; negative stress-test changes are retained in the CSV.
- Five repeated CIDs cover ten rows; cation/salt and metiram/zineb representation cases are distinguished from the other three structural conflicts.

This report does not claim that the full-D/exclude performance has been computed. Any separate fit is documented in its own sensitivity outputs.
