"""Verify nesting of the exported crossover test labels, without fitting.

The negative handling policy is an intentionally non-nested stress test;
it is not a feasible correction of three nested toxicity cutoffs. These
checks concern exported test cohorts; training labels are not exported
in the prediction file and are not claimed to have been checked here.
"""
from pathlib import Path
import hashlib
import json

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/label_audit"
SOURCE = OUT / "09_crossover_predictions.csv.gz"


def main():
    predictions = pd.read_csv(SOURCE)
    keys = ["Scenario", "Unresolved", "Split", "Index", "Threshold"]
    representations = sorted(predictions.Representation.unique())
    assert len(representations) == 5
    assert not predictions.duplicated(keys + ["Representation"]).any()
    grouped = predictions.groupby(keys, sort=True)
    assert grouped.size().eq(5).all(), "every label must be present for all five representations"
    for field in ["OriginalLabel", "CorrectedLabel", "CAS", "n_train", "n_test", "train_index_sha256", "test_index_sha256"]:
        assert grouped[field].nunique(dropna=False).eq(1).all(), f"representation-specific labels/membership: {field}"
    data = predictions[predictions.Representation == representations[0]]
    rows, examples = [], []
    for (scenario, handling, split), g in data.groupby(["Scenario", "Unresolved", "Split"], sort=True):
        assert set(g.Threshold) == {1, 11, 100}
        assert g.groupby("Index").size().eq(3).all()
        assert g.n_test.nunique() == 1 and g.n_train.nunique() == 1
        assert g.Index.nunique() == g.n_test.iloc[0]
        row = {"Scenario": scenario, "Unresolved": handling, "Split": split,
               "n_test": int(g.n_test.iloc[0]), "n_train": int(g.n_train.iloc[0]),
               "train_index_sha256": g.train_index_sha256.iloc[0],
               "test_index_sha256": g.test_index_sha256.iloc[0]}
        for field in ["OriginalLabel", "CorrectedLabel"]:
            matrix = g.pivot(index="Index", columns="Threshold", values=field).sort_index()
            assert matrix.isin([0, 1]).all().all()
            bad_1_11 = matrix[1] > matrix[11]
            bad_11_100 = matrix[11] > matrix[100]
            bad = bad_1_11 | bad_11_100
            row[f"{field}_nonnested"] = int(bad.sum())
            row[f"{field}_le1_exceeds_le11"] = int(bad_1_11.sum())
            row[f"{field}_le11_exceeds_le100"] = int(bad_11_100.sum())
            row[f"{field}_nonnested_indices"] = "|".join(map(str, matrix.index[bad]))
            for idx, label in matrix[bad].iterrows():
                identity = g[g.Index == idx].iloc[0]
                examples.append({"Scenario": scenario, "Unresolved": handling, "Split": split,
                                 "label_source": field, "Index": int(idx), "CAS": identity.CAS,
                                 "y_1": int(label[1]), "y_11": int(label[11]), "y_100": int(label[100])})
            if field == "OriginalLabel" or handling in ("keep", "exclude"):
                assert not bad.any(), (scenario, handling, split, field, matrix[bad])
        row["interpretation"] = ("non-nested extreme stress test; exclude from inference about feasible nested corrections"
                                 if handling == "negative" else "nested test labels")
        rows.append(row)
    counts = pd.DataFrame(rows)
    assert len(counts) == 45
    counts.to_csv(OUT / "LABEL_NESTING_CHECK.csv", index=False, lineterminator="\n")
    report = {"status": "passed_with_declared_non_nested_stress_tests",
              "scope": "exported test labels only; training labels are not in the prediction artifact",
              "input_file": str(SOURCE.relative_to(ROOT)).replace("\\", "/"),
              "input_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              "five_representation_label_and_membership_identity": True,
              "cohorts_checked": len(counts), "original_test_label_violations": int(counts.OriginalLabel_nonnested.sum()),
              "keep_exclude_test_label_violations": int(counts[counts.Unresolved != "negative"].CorrectedLabel_nonnested.sum()),
              "negative_policy_cohorts_with_violations": int((counts[counts.Unresolved == "negative"].CorrectedLabel_nonnested > 0).sum()),
              "counts": rows, "all_violating_test_identities": examples,
              "interpretation": "Violations under the negative policy are disclosed stress-test semantics, not evidence for a physically coherent corrected endpoint benchmark. Repeated compounds across splits/scenarios are not independent observations."}
    (OUT / "LABEL_NESTING_CHECK.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(counts[["Scenario", "Unresolved", "Split", "n_test", "OriginalLabel_nonnested", "CorrectedLabel_nonnested"]].to_string(index=False))


if __name__ == "__main__":
    main()
