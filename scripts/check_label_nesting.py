"""Reconstruct and check crossover labels directly from scientific inputs.

This check does not read a previous validation report or fit a model. It uses
the compound-level rule ledger, route-matched regulatory labels, and official
split files. A report is written only when --report-dir is requested.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
THRESHOLDS = (100, 11, 1)
SCENARIOS = ("primary", "conservative", "adult48", "regulatory", "combined")
POLICIES = ("keep", "exclude", "negative")
SPLITS = ("Random", "MaxMin", "Time")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def reconstruct_labels(root: Path):
    """Return original labels, pre-policy scenario labels and official indices."""
    raw = pd.read_csv(root / "data/raw/apistox.csv")
    require(len(raw) == 1035 and raw.SMILES.is_unique,
            "Expected 1,035 unique benchmark structures")
    tier = np.select(
        [raw.ppdb_level.eq(0), raw.ppdb_level.eq(1) & raw.label.eq(0),
         raw.ppdb_level.eq(1) & raw.label.eq(1), raw.ppdb_level.eq(2)],
        [0, 1, 2, 3], default=-1,
    )
    require((tier >= 0).all(), "Unmapped distributed label combination")
    original = np.array([(tier >= level).astype(np.int8) for level in (1, 2, 3)])
    index_by_smiles = dict(zip(raw.SMILES, raw.index))
    splits = {}
    for split in SPLITS:
        partitions = {}
        for partition, expected in (("train", 828), ("test", 207)):
            source = pd.read_csv(root / f"data/official_splits/{split.lower()}_{partition}.csv")
            mapped = source.SMILES.map(index_by_smiles)
            require(mapped.notna().all() and mapped.is_unique and len(mapped) == expected,
                    f"Invalid {split} {partition} membership")
            partitions[partition] = mapped.to_numpy(dtype=int)
        require(not set(partitions["train"]) & set(partitions["test"]),
                f"Overlapping {split} train/test identities")
        require(set(partitions["train"]) | set(partitions["test"]) == set(raw.index),
                f"Incomplete {split} membership")
        splits[split] = partitions

    ledger = pd.read_csv(root / "results/label_audit/02_compound_labels_by_rule.csv")
    ecotox_indices = set(raw.index[raw.source.eq("ECOTOX")])
    scenarios = {}
    for name, record_filter, rules in (
        ("primary", "all", ("D", "C", "D")),
        ("conservative", "all", ("C", "C", "C")),
        ("adult48", "adult48", ("D", "C", "D")),
    ):
        labels = original.copy()
        for axis, (threshold, rule) in enumerate(zip(THRESHOLDS, rules)):
            rows = ledger[(ledger["filter"] == record_filter)
                          & (ledger.routes == "eligible")
                          & (ledger.threshold == threshold) & (ledger.rule == rule)]
            require(rows.Index.is_unique and set(rows.Index) == ecotox_indices,
                    f"Incomplete or duplicate {name}/{threshold} rule rows")
            indices = rows.Index.to_numpy(dtype=int)
            require(np.array_equal(rows.original.to_numpy(), original[axis, indices]),
                    f"Misaligned original labels in {name}/{threshold}")
            require(rows.label.isin([-2, -1, 0, 1]).all(), "Unexpected rule outcome code")
            # -2 denotes no eligible record; policy handling treats it as unresolved.
            labels[axis, indices] = np.where(rows.label == -2, -1, rows.label)
        scenarios[name] = labels

    external = root / "output/external_data_audit"
    structures = pd.read_csv(external / "apistox_structures.csv")
    mapped = structures.SMILES.map(index_by_smiles)
    require(len(structures) == len(raw) and mapped.notna().all() and mapped.is_unique,
            "Regulatory structure map does not cover the benchmark")
    key_to_indices = {}
    for key, idx in zip(structures.inchikey, mapped.astype(int)):
        key_to_indices.setdefault(key, []).append(idx)
    concordance = pd.read_csv(external / "same_route_label_concordance.csv")
    concordance = concordance[concordance.source.isin(["OpenFoodTox", "PLOS2022"])]
    regulatory = np.full_like(original, -9)
    for key, rows in concordance.groupby("inchikey"):
        require(key in key_to_indices and len(key_to_indices[key]) == 1,
                f"Ambiguous or missing regulatory identity: {key}")
        idx = key_to_indices[key][0]
        for axis, threshold in enumerate(THRESHOLDS):
            require((rows[f"apistox_y_{threshold}"] == original[axis, idx]).all(),
                    "Regulatory rows do not match distributed labels")
            values = set(rows[f"external_y_{threshold}"].astype(int))
            require(values <= {-1, 0, 1}, "Unexpected regulatory label code")
            regulatory[axis, idx] = next(iter(values)) if len(values) == 1 else -1
    scenarios["regulatory"] = np.where(regulatory == -9, original, regulatory)
    combined = scenarios["regulatory"].copy()
    is_ecotox = raw.source.eq("ECOTOX").to_numpy()
    for axis in range(3):
        has_regulatory = regulatory[axis] != -9
        only_ecotox = is_ecotox & ~has_regulatory
        both = is_ecotox & has_regulatory
        combined[axis, only_ecotox] = scenarios["primary"][axis, only_ecotox]
        combined[axis, both] = np.where(
            scenarios["primary"][axis, both] == regulatory[axis, both],
            regulatory[axis, both], -1,
        )
    scenarios["combined"] = combined
    return original, scenarios, splits


def apply_policy(labels, original, policy):
    retained = (labels >= 0).all(axis=0) if policy == "exclude" else np.ones(labels.shape[1], bool)
    resolved = np.where(labels < 0, original, labels)
    if policy == "negative":
        resolved[0] = np.where(labels[0] < 0, 0, labels[0])
    return resolved, retained


def violations(labels):
    return (labels[2] > labels[1]) | (labels[1] > labels[0])


def check_label_nesting(root: Path = ROOT, check_predictions: bool = True):
    """Validate 45 paired cohorts and return a summary plus 90 partition rows.

Prediction validation compares every representation's exported labels and test
membership with the reconstructed labels, without using old pass/fail outputs.
The negative policy is a declared nonnested stress test, not a failure.
"""
    original, scenarios, splits = reconstruct_labels(root)
    require(not violations(original).any(), "Distributed labels are not nested")
    predictions = None
    if check_predictions:
        predictions = pd.read_csv(root / "results/label_audit/09_crossover_predictions.csv.gz")
        require(set(predictions.Scenario) == set(SCENARIOS), "Prediction scenario coverage differs")
        require(set(predictions.Unresolved) == set(POLICIES), "Prediction policy coverage differs")
        require(set(predictions.Split) == set(SPLITS), "Prediction split coverage differs")
        require(set(predictions.Representation) == {"ECFP", "Avalon", "MACCS", "WL-HI", "Descriptors12"},
                "Prediction representation coverage differs")
    records = []
    checked_predictions = 0
    for scenario in SCENARIOS:
        for policy in POLICIES:
            labels, retained = apply_policy(scenarios[scenario], original, policy)
            for split, partitions in splits.items():
                members = {part: ix[retained[ix]] for part, ix in partitions.items()}
                for partition, indices in members.items():
                    bad = violations(labels[:, indices])
                    require(policy == "negative" or not bad.any(),
                            f"Nesting violation in {scenario}/{policy}/{split}/{partition}")
                    records.append({
                        "Scenario": scenario, "Unresolved": policy, "Split": split,
                        "Partition": partition, "n": len(indices),
                        "OriginalLabel_nonnested": int(violations(original[:, indices]).sum()),
                        "CorrectedLabel_nonnested": int(bad.sum()),
                        "CorrectedLabel_nonnested_indices": ";".join(map(str, indices[bad])),
                    })
                if predictions is not None:
                    cohort = predictions[(predictions.Scenario == scenario)
                                         & (predictions.Unresolved == policy)
                                         & (predictions.Split == split)]
                    indices = members["test"]
                    require(len(cohort) == 15 * len(indices), "Incomplete crossover prediction cohort")
                    for (representation, threshold), group in cohort.groupby(["Representation", "Threshold"]):
                        require(threshold in THRESHOLDS, "Unexpected prediction threshold")
                        require(group.Index.is_unique and set(group.Index) == set(indices),
                                f"Prediction membership differs: {scenario}/{policy}/{split}/{representation}/{threshold}")
                        require((group.n_train == len(members["train"])).all()
                                and (group.n_test == len(indices)).all(), "Prediction cohort sizes differ")
                        axis = THRESHOLDS.index(threshold)
                        ordered = group.set_index("Index").loc[indices]
                        require(np.array_equal(ordered.OriginalLabel, original[axis, indices])
                                and np.array_equal(ordered.CorrectedLabel, labels[axis, indices]),
                                f"Exported labels differ: {scenario}/{policy}/{split}/{representation}/{threshold}")
                        checked_predictions += len(group)
    table = pd.DataFrame(records)
    summary = {
        "status": "passed_with_declared_nonnested_stress_tests",
        "method": "Direct reconstruction from compound rule ledger, regulatory concordance, and official splits",
        "cohorts": 45, "partition_checks": len(table),
        "exported_prediction_rows_checked": checked_predictions,
        "original_label_violations": int(table.OriginalLabel_nonnested.sum()),
    }
    for partition in ("train", "test"):
        rows = table[table.Partition == partition]
        negative = rows[rows.Unresolved == "negative"]
        summary[f"keep_exclude_{partition}_violations"] = int(
            rows[rows.Unresolved.isin(["keep", "exclude"])].CorrectedLabel_nonnested.sum())
        summary[f"negative_policy_{partition}_cohorts_with_violations"] = int(
            (negative.CorrectedLabel_nonnested > 0).sum())
        summary[f"negative_policy_{partition}_violations_min"] = int(negative.CorrectedLabel_nonnested.min())
        summary[f"negative_policy_{partition}_violations_max"] = int(negative.CorrectedLabel_nonnested.max())
    return summary, table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT, help="Research checkout to inspect")
    parser.add_argument("--report-dir", type=Path, help="Optional directory for JSON and CSV reports")
    args = parser.parse_args()
    summary, rows = check_label_nesting(args.root.resolve())
    if args.report_dir is not None:
        args.report_dir.mkdir(parents=True, exist_ok=True)
        (args.report_dir / "label_nesting_summary.json").write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        rows.to_csv(args.report_dir / "label_nesting_partitions.csv", index=False)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
