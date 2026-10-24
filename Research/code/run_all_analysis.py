"""Regenerate the numerical analyses from the supplied scientific inputs.

This entry point regenerates all numerical analysis outputs from the frozen raw
CSV, six official split files, and cached source records. Manuscript compilation
and figure rendering are separate steps. Label checks are computed directly;
pre-existing validation reports are not required.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
RESULTS = ROOT / "results"
SCRIPTS = [
    "00_validate_and_describe.py",
    "01_primary_endpoint_models.py",
    "02_metadata_decomposition.py",
    "03_matched_learning_curves.py",
    "04_conditional_neighborhood.py",
    "05_scaffold_sensitivity.py",
    "06_novelty_calibration.py",
    "07_sparse_descriptors_supplementary.py",
    "08_sensitivity_diagnostics.py",
    "09_export_submission_tables.py",
    "10_tier_boundary_diagnostics.py",
    "11_qualifier_propagation_audit.py",
    "12_label_rule_audit.py",
    "13_corrected_label_crossover.py",
    "14_record_case_audit.py",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clean", action="store_true", help="Remove results before running.")
    parser.add_argument("--validation-report-dir", type=Path,
                        help="Optional output directory for label-check reports.")
    args = parser.parse_args()
    if args.clean and RESULTS.exists():
        if RESULTS.resolve() != ROOT.resolve() / "results":
            raise ValueError("Refusing to clean results outside this checkout")
        shutil.rmtree(RESULTS)
    RESULTS.mkdir(parents=True, exist_ok=True)
    for script in SCRIPTS:
        path = CODE / script
        print(f"\n=== {script} ===", flush=True)
        subprocess.run([sys.executable, str(path)], cwd=ROOT, check=True)
    subprocess.run([sys.executable, str(ROOT / "scripts" / "prepare_cache.py")],
                   cwd=ROOT, check=True)
    subprocess.run([sys.executable, str(ROOT / "scripts" / "describe_cohort_selection.py")],
                   cwd=ROOT, check=True)
    check_command = [sys.executable, str(ROOT / "scripts" / "check_label_nesting.py")]
    if args.validation_report_dir is not None:
        check_command += ["--report-dir", str(args.validation_report_dir.resolve())]
    subprocess.run(check_command, cwd=ROOT, check=True)
    print("\nAll numerical analyses and direct label checks completed.")


if __name__ == "__main__":
    main()
