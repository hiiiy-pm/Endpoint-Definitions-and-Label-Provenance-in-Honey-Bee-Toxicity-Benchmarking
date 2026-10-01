"""Run the final numerical analyses for the English submission package.

This entry point regenerates all numerical analysis outputs from the frozen raw
CSV and the six official split files. Manuscript compilation and figure rendering
are handled separately by build.sh and scripts/.
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
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clean", action="store_true", help="Remove results before running.")
    args = parser.parse_args()
    if args.clean and RESULTS.exists():
        shutil.rmtree(RESULTS)
    RESULTS.mkdir(parents=True, exist_ok=True)
    for script in SCRIPTS:
        path = CODE / script
        print(f"\n=== {script} ===", flush=True)
        subprocess.run([sys.executable, str(path)], cwd=ROOT, check=True)
    print("\nAll final numerical analyses completed.")


if __name__ == "__main__":
    main()
