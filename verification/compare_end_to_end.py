"""Compare a clean numerical rerun with the working analysis, preserving evidence."""
from __future__ import annotations

import argparse
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ATOL = 1e-12
OPTIONAL = {"results/record_audit/INDEPENDENT_VALIDATION.json",
            "results/record_audit/UPSTREAM_PIN_CHECK.json"}
COMPARISON_OUTPUTS = {"end_to_end_comparison.csv", "END_TO_END_VERIFICATION.json",
                      "end_to_end_run.txt", "end_to_end_nesting_run.txt", "end_to_end_structure_run.txt"}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def array_difference(a, b):
    if a.shape != b.shape:
        raise AssertionError(f"Array shapes differ: {a.shape}, {b.shape}")
    if np.issubdtype(a.dtype, np.number) and np.issubdtype(b.dtype, np.number):
        if not np.array_equal(np.isnan(a), np.isnan(b)):
            raise AssertionError("Numeric missing-value patterns differ")
        if not np.allclose(a, b, atol=ATOL, rtol=0, equal_nan=True):
            raise AssertionError(f"Numeric values differ beyond atol={ATOL}")
        finite = np.isfinite(a) & np.isfinite(b)
        return float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else 0.
    if not np.array_equal(a, b):
        raise AssertionError("Nonnumeric array values differ")
    return 0.


def csv_difference(a_path, b_path):
    a = pd.read_csv(a_path, float_precision="round_trip")
    b = pd.read_csv(b_path, float_precision="round_trip")
    if a.shape != b.shape or list(a.columns) != list(b.columns):
        raise AssertionError(f"CSV shape/columns differ: {a.shape}, {b.shape}")
    if not np.array_equal(a.isna().to_numpy(), b.isna().to_numpy()):
        raise AssertionError("CSV missing-value patterns differ")
    difference = 0.
    for column in a:
        if pd.api.types.is_numeric_dtype(a[column]) and pd.api.types.is_numeric_dtype(b[column]):
            difference = max(difference, array_difference(a[column].to_numpy(), b[column].to_numpy()))
        elif not a[column].fillna("<NA>").astype(str).equals(b[column].fillna("<NA>").astype(str)):
            raise AssertionError(f"CSV text differs: {column}")
    return difference, str(a.shape)


def json_difference(a, b, path=""):
    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            raise AssertionError(f"JSON keys differ at {path}")
        return max([json_difference(a[k], b[k], f"{path}/{k}") for k in a] or [0.])
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            raise AssertionError(f"JSON list lengths differ at {path}")
        return max([json_difference(x, y, f"{path}/{i}") for i, (x, y) in enumerate(zip(a, b))] or [0.])
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not np.isclose(a, b, atol=ATOL, rtol=0, equal_nan=True):
            raise AssertionError(f"JSON numeric values differ at {path}: {a}, {b}")
        return abs(float(a) - float(b)) if np.isfinite(a) and np.isfinite(b) else 0.
    if a != b:
        raise AssertionError(f"JSON value differs at {path}: {a!r}, {b!r}")
    return 0.


def normalize_metadata(path, value):
    value = copy.deepcopy(value)
    excluded = []
    def remove(keys):
        for key in keys:
            if key in value:
                value.pop(key)
                excluded.append(key)
    if path == "results/cache/CACHE_METADATA.json":
        # The legacy record has no historical parameters/environment. Those fields
        # are validated against the current raw-SMILES attestation instead.
        value = {k: value[k] for k in ("seed", "n_molecules", "source_sha256")}
        excluded = ["legacy schema versus freshly generated schema-v2 provenance"]
    elif path == "results/cache/CACHE_VALIDATION.json":
        remove(["validated_at_utc", "legacy_generation_environment", "legacy_metadata_sha256"])
        value["artifact_sha256"] = sorted(value["artifact_sha256"])
        excluded.append("artifact byte hashes (all cache array/CSV contents compared separately)")
    elif path == "results/final/00_data_validation.json":
        remove(["cache_validation"])
    elif path == "results/final/EXPORT_PROVENANCE.json":
        value["input_sha256"] = sorted(value["input_sha256"])
        for output in value["outputs"].values():
            output.pop("sha256", None)
        excluded = ["input/output byte hashes (scientific values compared separately)"]
    elif path == "results/reproducibility/revision_verification.json":
        value["cache"], ignored = normalize_metadata("results/cache/CACHE_VALIDATION.json", value["cache"])
        value["negative_paths"] = [case for case in value["negative_paths"] if case["case"] != "legacy_without_validation"]
        excluded = ["legacy-only negative test absent for schema-v2 cache", *ignored]
    elif path == "output/descriptor_audit/protocol_and_verification.json":
        remove(["elapsed_seconds", "original_input_sha256", "script_sha256", "output_sha256"])
        value["versions"].pop("python")
        excluded.append("historical Python 3.12.14 versus rerun Python 3.12.13; package versions unchanged")
    elif path == "results/label_audit/LABEL_NESTING_CHECK.json":
        remove(["input_sha256"])
    elif path == "results/record_audit/SELECTION_STRUCTURE_METADATA.json":
        value["input_sha256"] = sorted(value["input_sha256"])
        excluded.append("input byte hashes, including legacy versus fresh cache attestation (input contents compared separately)")
    return value, excluded


def files(root):
    output = {}
    for folder in ("results", "source_data", "output/descriptor_audit"):
        for path in (root / folder).rglob("*"):
            if path.is_file() and path.name not in COMPARISON_OUTPUTS:
                output[path.relative_to(root).as_posix()] = path
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("rerun", type=Path)
    args = parser.parse_args()
    rerun = args.rerun.resolve()
    execution = json.loads((rerun / "e2e_execution.json").read_text(encoding="utf-8-sig"))
    supplemental_execution = json.loads((rerun / "nesting_execution.json").read_text(encoding="utf-8-sig"))
    structure_execution = json.loads((rerun / "structure_execution.json").read_text(encoding="utf-8-sig"))
    if any(record["exit_code"] != 0 for record in (execution, supplemental_execution, structure_execution)):
        raise ValueError("The numerical rerun did not exit successfully")
    current, fresh = files(ROOT), files(rerun)
    regenerated_source = set(json.loads((rerun / "results/final/EXPORT_PROVENANCE.json").read_text())["outputs"])
    regenerated_source.update("source_data/" + name for name in (
        "figure_tier_pair_auroc.csv", "figure_tier1_exclusion.csv", "figure_concordance_by_tier.csv", "figure_qualifier_propagation.csv"))
    rows = []
    for relative in sorted(current.keys() | fresh.keys()):
        a, b = current.get(relative), fresh.get(relative)
        scope = "regenerated"
        if relative in OPTIONAL:
            scope = "optional_separate_verification"
        elif relative.startswith("source_data/") and relative not in regenerated_source:
            scope = "retained_external_or_separate_figure_output"
        row = {"path": relative, "scope": scope, "status": "", "bytes_identical": False,
               "max_absolute_difference": 0., "shape": "", "classification_note": "",
               "working_sha256": sha(a) if a else "", "rerun_sha256": sha(b) if b else ""}
        if a is None or b is None:
            row["status"] = "not_regenerated" if scope != "regenerated" else "FAIL"
            row["classification_note"] = "Absent from " + ("working repository" if a is None else "clean rerun")
            rows.append(row)
            continue
        row["bytes_identical"] = row["working_sha256"] == row["rerun_sha256"]
        try:
            if relative.endswith((".csv", ".csv.gz")):
                row["max_absolute_difference"], row["shape"] = csv_difference(a, b)
            elif a.suffix == ".npy":
                aa, bb = np.load(a, allow_pickle=False), np.load(b, allow_pickle=False)
                row["max_absolute_difference"] = array_difference(aa, bb)
                row["shape"] = str(aa.shape)
            elif a.suffix == ".npz":
                with np.load(a, allow_pickle=False) as aa, np.load(b, allow_pickle=False) as bb:
                    if aa.files != bb.files:
                        raise AssertionError("NPZ keys differ")
                    row["max_absolute_difference"] = max(array_difference(aa[key], bb[key]) for key in aa.files)
                    row["shape"] = f"{len(aa.files)} arrays"
            elif a.suffix == ".json":
                aa, exclusions = normalize_metadata(relative, json.loads(a.read_text(encoding="utf-8-sig")))
                bb, _ = normalize_metadata(relative, json.loads(b.read_text(encoding="utf-8-sig")))
                row["max_absolute_difference"] = json_difference(aa, bb)
                row["classification_note"] = "; ".join(exclusions)
            else:
                if a.read_text(encoding="utf-8-sig").replace("\r\n", "\n") != b.read_text(encoding="utf-8-sig").replace("\r\n", "\n"):
                    raise AssertionError("Text differs after newline normalization")
            row["status"] = ("byte_identical" if row["bytes_identical"] else
                             "metadata_with_numeric_roundoff" if row["classification_note"] and row["max_absolute_difference"] else
                             "metadata_only" if row["classification_note"] else
                             "numeric_tolerance" if row["max_absolute_difference"] else "value_identical")
        except (AssertionError, ValueError) as error:
            row["status"] = "FAIL"
            row["classification_note"] = str(error)
        rows.append(row)
    output = ROOT / "results/reproducibility"
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output / "end_to_end_comparison.csv", index=False)
    failures = [row for row in rows if row["status"] == "FAIL"]
    generated = [row for row in rows if row["scope"] == "regenerated"]
    code_inputs = [*sorted((rerun / "code").glob("*.py")), rerun / "verification/verify_revision.py",
                   rerun / "verification/check_crossover_nesting.py",
                   rerun / "verification/describe_cohort_selection.py",
                   rerun / "requirements-analysis.txt"]
    code_changes = [path.relative_to(rerun).as_posix() for path in code_inputs
                    if sha(path) != sha(ROOT / path.relative_to(rerun))]
    report = {"status": "PASS" if not failures and not code_changes else "FAIL", "execution": execution,
              "supplemental_nesting_execution": supplemental_execution,
              "supplemental_structure_execution": structure_execution,
              "working_repository": str(ROOT), "isolated_rerun": str(rerun),
              "scripts_completed": [f"{i:02d}" for i in range(15)], "revision_verifier_completed": True,
              "crossover_nesting_verifier_completed": True,
              "cohort_structure_verifier_completed": True,
              "cohort_structure_independent_raw_smiles_maxima_checked": 3105,
              "entrypoint_update_after_numerical_run": "check_crossover_nesting.py and describe_cohort_selection.py were appended to run_all after the 00-14 run and executed separately on its fresh outputs; numerical scripts were unchanged",
              "initial_numerical_entrypoint_sha256": sha(rerun / "run_all_executed_initial.py"),
              "environment": json.loads((rerun / "results/reproducibility/revision_verification.json").read_text())["environment"],
              "comparison_policy": {"numeric_absolute_tolerance": ATOL, "numeric_relative_tolerance": 0,
                                    "CSV": "same shapes, columns, row order, text and missing-value pattern",
                                    "arrays": "same shapes/keys and missing-value patterns; all stored bootstrap arrays included"},
              "compared_file_count": len(rows), "regenerated_file_count": len(generated),
              "status_counts": dict(Counter(row["status"] for row in rows)),
              "regenerated_status_counts": dict(Counter(row["status"] for row in generated)),
              "max_absolute_scientific_difference": max(row["max_absolute_difference"] for row in generated),
              "failures": failures, "analysis_code_changed_since_copy": code_changes,
              "analysis_code_sha256": {path.relative_to(rerun).as_posix(): sha(path) for path in code_inputs},
              "separate_not_regenerated": [row["path"] for row in rows if row["scope"] != "regenerated"],
              "log_sha256": sha(rerun / "e2e_run.log"),
              "nesting_log_sha256": sha(rerun / "nesting_run.log"),
              "structure_log_sha256": sha(rerun / "structure_run.log"),
              "historical_descriptor_python": json.loads((ROOT / "output/descriptor_audit/protocol_and_verification.json").read_text())["versions"]["python"],
              "current_descriptor_python": json.loads((rerun / "output/descriptor_audit/protocol_and_verification.json").read_text())["versions"]["python"],
              "limits": ["The optional upstream network pin check and separate record-workbook validation are not rerun by the offline numerical entry point.",
                         "Separate figure/external source_data artifacts retained during copying are not counted as regenerated.",
                         "Manuscript compilation, figure rendering and external study acquisition are outside this numerical rerun."]}
    (output / "END_TO_END_VERIFICATION.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    shutil.copyfile(rerun / "e2e_run.log", output / "end_to_end_run.txt")
    shutil.copyfile(rerun / "nesting_run.log", output / "end_to_end_nesting_run.txt")
    shutil.copyfile(rerun / "structure_run.log", output / "end_to_end_structure_run.txt")
    print(json.dumps({key: report[key] for key in ("status", "execution", "status_counts", "regenerated_status_counts", "max_absolute_scientific_difference", "analysis_code_changed_since_copy", "failures")}, indent=2))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
