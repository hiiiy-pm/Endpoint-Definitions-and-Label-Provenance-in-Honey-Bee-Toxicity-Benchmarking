# Clean end-to-end numerical rerun

The current revision was rerun from raw SMILES and frozen source exports in a new
isolated repository copy at `D:\ApisTox_revision_e2e_20261002`. The directory name
retains the revision date; the actual run was **3 October 2026, 00:00:38–00:04:43
Asia/Shanghai**, taking 244.96 seconds. All scripts **00–14** and
`verification/verify_revision.py` completed with exit code zero.

The exported crossover-label nesting check was then appended to the entry point
and executed separately against the fresh predictions at 00:06:34–00:06:35
(1.20 seconds, exit zero). No numerical script changed between the clean run and
this additional check. The original entry point's hash and this two-step execution
history are retained in the machine-readable report. The subsequent descriptive
cohort-structure check was also appended to the entry point and executed against
the fresh cache and record cohorts at 00:11:59–00:12:00 (1.08 seconds, exit zero).
It independently checked 3,105 raw-SMILES similarity maxima and regenerated its
two CSVs and metadata. Future invocations run all three verification programs
automatically. Neither additional check required repeating model training.

## Environment and command

The run used the existing pinned interpreter
`C:\Users\hiy\.venvs\apistox\Scripts\python.exe`: Python 3.12.13, NumPy 2.5.3,
pandas 3.0.1, SciPy 1.18.1, scikit-learn 1.8.0 and RDKit 2025.9.4.
The scientific dependencies match `requirements-analysis.txt`.

```shell
python -u code/run_all_analysis.py --clean
```

The isolated copy retained the external source inputs but excluded `.git`,
`__pycache__` and ignored `verification/results`. `--clean` removed only that
copy's `results` directory after checking its resolved location. All four molecular
kernels, 12 descriptors and scaffold identifiers were regenerated from the raw
SMILES. The working repository's legacy cache was not replaced.

The older `output/descriptor_audit/protocol_and_verification.json` records Python
3.12.14 for that historical selector run. This rerun used Python 3.12.13 and
reproduced its scientific outputs. Both environment records are retained; the
older record is not presented as evidence of this run.

## File comparison

`verification/compare_end_to_end.py` compared 144 files across `results`,
`source_data` and `output/descriptor_audit`. Of these, **140 were regenerated** by
the numerical pipeline and its verification steps:

| Comparison of regenerated files | Files |
|---|---:|
| Byte-identical | 95 |
| Different bytes, identical scientific values/text | 24 |
| Floating-point differences within the declared tolerance | 14 |
| Only declared provenance/environment differences | 6 |
| Declared provenance differences plus numerical roundoff in a verification summary | 1 |

Every compared scientific CSV has the same shape, column names, row order,
text values and missing-value pattern. Numeric comparison uses absolute tolerance
**1e-12 and zero relative tolerance**. The maximum difference is
**4.159450561758149e-13**, in model-decomposition decision scores. NPY/NPZ shapes,
array keys, values and missingness were checked directly, including all saved
bootstrap and permutation draws. These arrays were numerically identical.

The regenerated primary bootstrap intervals, crossover AUROCs/gaps/intervals,
label-rule and record-case tables, and the scientific tables exported by script 09
agree with the working results under this policy. Nested descriptor selection
again fitted 54 selectors and recorded 34,668 fit events, with zero forbidden-row
events. Crossover checks independently reconstructed the retained `exclude`
cohorts, refitted both training-label models in 225 cells (450 fits), and verified
all 2,700 AUROC and 2,700 gap rows from the 133,545 saved predictions.

The two cohort-selection structure tables were byte-identical. The metadata
differs only in byte hashes of upstream artifacts, including the legacy versus
new cache attestation. Direct RDKit similarity maxima agreed with the float32
cache within 2.95e-8; that representation-precision comparison is separate from
the exact agreement of the working and clean-rerun structure tables. These maxima
describe similarity to the selected reference sets and do not establish sample
representativeness.

The label-nesting verifier checked 45 test cohorts. Original labels and corrected
`keep`/`exclude` test labels are nested; all 15 `negative` stress-test cohorts
contain disclosed nesting violations. This extreme policy does not define a
physically coherent corrected cumulative-endpoint benchmark. The nesting artifact
itself checks exported test labels, not unexported training labels.

Differences in cache-generation schema, timestamps, input/output byte hashes,
line endings and historical environment metadata are classified separately.
Scientific array and CSV comparison is not bypassed by this classification, and
the cache validator's identity, hash, shape and environment checks remain strict.

## Evidence and scope

- `results/reproducibility/END_TO_END_VERIFICATION.json`: run times, environment,
  code hashes, comparison policy, counts and zero failures.
- `results/reproducibility/end_to_end_comparison.csv`: per-file hashes,
  regeneration scope, classification and maximum numeric differences.
- `results/reproducibility/end_to_end_run.txt`: complete clean numerical-run log.
- `results/reproducibility/end_to_end_nesting_run.txt`: additional nesting-check log.
- `results/reproducibility/end_to_end_structure_run.txt`: additional cohort-structure-check log.
- `results/reproducibility/revision_verification.json`: working-cache verification;
  the isolated copy contains the corresponding freshly generated-cache report.

To repeat the comparison against this retained isolated run:

```shell
python verification/compare_end_to_end.py D:\ApisTox_revision_e2e_20261002
```

Four files are explicitly outside the regenerated-file count. The separate
`record_audit/UPSTREAM_PIN_CHECK.json` network check and
`record_audit/INDEPENDENT_VALIDATION.json` original-workbook validation are not
produced by the offline numerical entry point. The external-reassessment and
label-rule-comparison figure CSVs under `source_data` were retained during copying
and were byte-identical, but are not claimed to have been regenerated by scripts
00–14. Manuscript compilation, figure rendering, and acquiring new external
studies are outside this numerical rerun.
