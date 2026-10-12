"""Assemble the revision statistics report from completed verification outputs."""
from pathlib import Path
import hashlib
import json

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]


def main():
    names = ["STATISTICS_VERIFICATION.json", "UNIFORM_D_EXCLUDE_METADATA.json",
             "CURRENT_IDENTIFIER_FLAGGED_EXCLUSION_METADATA.json", "DERIVED_OUTPUT_VERIFICATION.json",
             "figure_qa/FIGURE_NUMERIC_METADATA.json", "figure_qa/Fig7_label_crossover.alignment.json",
             "figure_qa/Fig7_label_crossover.collision.json"]
    payload = {name: json.loads((OUT / name).read_text(encoding="utf-8")) for name in names}
    assert payload["STATISTICS_VERIFICATION.json"]["Status"] == "passed"
    assert payload["DERIVED_OUTPUT_VERIFICATION.json"]["Status"] == "passed"
    report_dir = ROOT.parent / "revision_audit"
    report_dir.mkdir(exist_ok=True)
    (report_dir / "statistics_report.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    report = """# Revision statistics and crossover figure

All requested statistics work is complete. Frozen data, representations and result files were read without alteration. New model outputs are confined to `research/supplementary_audits/statistics/`.

## Label and source-log replay

- All-D/keep and D/C/D/keep have identical three full-cohort label vectors. The only unresolved-rule differences are indices 861 and 953 at 11 µg/bee; keep resolves them identically. Label-vector byte hashes are retained.
- Exclude retains 913 versus 915 identities. The two additions are training members in Random and MaxMin and test members in Time.
- The 745-entry CAS-level uncertainty log contains 536 observations in excluded route groups and 209 in retained groups; 119 of the latter remain in final ECOTOX determining groups.
- A Benomyl operator replay changes provisional rule determinations but leaves all five reported scenarios' effective keep/exclude inputs unchanged. The archive is not edited and the primary-paper adjudication remains in the separate source audit.
- Five repeated CID groups contain ten rows. The inventory distinguishes the chlormequat cation/chloride case and the metiram/zineb representation case from three other structure conflicts.

## Completed model sensitivities

Uniform-D/exclude was fitted separately for ECFP in all three official splits, with all four training/evaluation label cells and both endpoint contrasts. The full original and mixed-D/C/D exclude fits passed 27 saved-score reproduction checks, with maximum absolute discrepancy 2.22e-16. All six Time fits with unchanged training inputs were reused. Both severe-minus-broad contrasts have the same full-revision change within each split: Random +0.011472 (95% CI -0.069860 to 0.098371), MaxMin +0.008936 (-0.066118 to 0.090622), Time -0.008346 (-0.098811 to 0.085287). These use 2,000 shared scaffold draws and condition on selected cohorts and fitted models.

The predefined current-identifier screening set excludes 33 flags, comprising six parent/fragment cases and 27 remaining differences requiring adjudication. These are not 33 confirmed historical errors. Its frozen CSV hash was checked before model fitting. Retained training/test counts are 804/198, 811/191 and 808/194 for Random, MaxMin and Time. Nine original fits reproduced saved scores to 2.22e-16, then nine reduced fits were run. Original/all, frozen/reduced and refit/reduced results are all retained. Refit AUROCs at 100/11/1 are 0.715833/0.855508/0.921028, 0.586522/0.798537/0.795807 and 0.673675/0.766419/0.786168. The broad endpoint remains lowest in each split. Its two MaxMin contrasts have simultaneous 95% CIs above zero; all paired intervals for total changes in the gaps include zero. This analysis uses 5,000 shared scaffold draws and two-contrast simultaneous intervals, distinct from the 2,000-draw label crossover.

An independent standard-library verifier recomputed 99 AUROCs and 78 gaps by enumerating all positive-negative pairs, without sklearn, and checked 183 membership hashes/prediction orders against the official splits and declared filters. All 177 numeric checks passed; maximum discrepancy was 2.22e-16.

## Main crossover figure

`Fig7_label_crossover.pdf` and `.png` contain all 24 primary/keep/ECFP four-cell gap estimates and all six full-revision paired effects. Source data are in `source_data/figure_label_crossover.csv`. Thirty independent numeric checks, including regeneration of the original bootstrap intervals, passed with maximum discrepancy 9.7e-17. Rendered alignment, PDF font-size and collision checks passed; the minimum glyph size is 7 pt, with zero final collision failures or warnings. Per-panel visual QA is recorded in `figure_qa/RENDERED_QA.md`.

## Reproduction

Run from the research directory using the pinned project Python environment:

```text
python supplementary_audits/statistics/verify_audit_statistics.py
python supplementary_audits/statistics/run_uniform_d_exclude.py
python supplementary_audits/statistics/run_current_identifier_flagged_exclusion.py
python supplementary_audits/statistics/verify_derived_outputs.py
python supplementary_audits/statistics/make_statistics_si.py
python scripts/make_crossover_figure.py
python supplementary_audits/statistics/finalize_statistics_report.py
```

The two independent verification scripts require only the Python standard library. The model and figure scripts use the pinned environment. Full figure QA commands and their output are retained under `figure_qa/`; the SI identifier-screen insert is generated directly from the saved CSV tables.
"""
    (report_dir / "statistics_report.md").write_text(report, encoding="utf-8")
    manifest = {str(p.relative_to(OUT)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(OUT.rglob("*")) if p.is_file() and "__pycache__" not in p.parts and p.name != "OUTPUT_SHA256.json"}
    (OUT / "OUTPUT_SHA256.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(report_dir / "statistics_report.md")


if __name__ == "__main__":
    main()
