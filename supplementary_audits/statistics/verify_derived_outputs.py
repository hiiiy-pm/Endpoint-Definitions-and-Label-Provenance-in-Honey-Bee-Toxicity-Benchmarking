"""Recompute newly fitted AUROCs and gaps by explicit pair counting, without sklearn."""
from __future__ import annotations
import csv
import hashlib
import json
import struct
from collections import defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]


def read(name):
    with (OUT / name).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def auc(labels, scores):
    pos = [s for y, s in zip(labels, scores) if y == 1]
    neg = [s for y, s in zip(labels, scores) if y == 0]
    assert pos and neg
    return sum((p > n) + .5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))


def main():
    results = []
    benchmark = read(ROOT / "data/raw/apistox.csv")
    sm2index = {r["SMILES"]: i for i, r in enumerate(benchmark)}
    original = {}
    for t, cutoff in ((100, 1), (11, 2), (1, 3)):
        original[t] = [int((0 if int(r["ppdb_level"]) == 0 else
                            3 if int(r["ppdb_level"]) == 2 else 1 + int(r["label"])) >= cutoff) for r in benchmark]
    official = {}
    for split in ("Random", "MaxMin", "Time"):
        for part in ("train", "test"):
            official[split, part] = [sm2index[r["SMILES"]] for r in read(ROOT / f"data/official_splits/{split.lower()}_{part}.csv")]
        assert not set(official[split, "train"]) & set(official[split, "test"])
        assert set(official[split, "train"]) | set(official[split, "test"]) == set(range(1035))
    audit = read(ROOT / "results/label_audit/02_compound_labels_by_rule.csv")
    excluded_uniform = {}
    for scenario in ("mixed_DCD", "uniform_DDD"):
        excluded_uniform[scenario] = {int(r["Index"]) for r in audit if r["filter"] == "all" and r["routes"] == "eligible"
                                     and r["rule"] == ("C" if scenario == "mixed_DCD" and r["threshold"] == "11" else "D")
                                     and int(r["label"]) < 0}
    def index_hash(indices):
        return hashlib.sha256(b"".join(struct.pack("<q", i) for i in indices)).hexdigest()
    membership_checks = 0
    predictions = defaultdict(list)
    for row in read("uniform_d_exclude_predictions.csv"):
        predictions[row["Scenario"], row["Split"], row["Threshold"]].append(row)
    estimates = {}
    for row in read("uniform_d_exclude_auroc.csv"):
        key = row["Scenario"], row["Split"], row["Threshold"]
        records = predictions[key]
        expected_test = [i for i in official[row["Split"], "test"] if i not in excluded_uniform[row["Scenario"]]]
        expected_train = [i for i in official[row["Split"], "train"] if i not in excluded_uniform[row["Scenario"]]]
        assert [int(r["Index"]) for r in records] == expected_test
        assert int(row["n_train"]) == len(expected_train) and int(row["n_test"]) == len(expected_test)
        assert row["train_index_sha256"] == index_hash(expected_train)
        assert row["test_index_sha256"] == index_hash(expected_test)
        membership_checks += 2
        yc = "OriginalLabel" if row["EvalLabels"] == "original" else "ProvisionalLabel"
        sc = "OriginalTrainScore" if row["TrainLabels"] == "original" else "ProvisionalTrainScore"
        value = auc([int(r[yc]) for r in records], [float(r[sc]) for r in records])
        error = abs(value - float(row["AUROC"]))
        assert error < 1e-12
        estimates[*key, row["TrainLabels"], row["EvalLabels"]] = value
        results.append({"Analysis": "uniform_D_exclude", "Quantity": "AUROC", "MaxAbsDifference": error})
    for row in read("uniform_d_exclude_gap.csv"):
        s, split, t = row["Scenario"], row["Split"], row["Severe"]
        a, b = row["TrainLabels"], row["EvalLabels"]
        if a.startswith("change:"):
            value = ((estimates[s, split, t, "provisional", "provisional"] - estimates[s, split, "100", "provisional", "provisional"])
                     - (estimates[s, split, t, "original", "original"] - estimates[s, split, "100", "original", "original"]))
        else:
            value = estimates[s, split, t, a, b] - estimates[s, split, "100", a, b]
        error = abs(value - float(row["Gap"]))
        assert error < 1e-12
        results.append({"Analysis": "uniform_D_exclude", "Quantity": "Gap", "MaxAbsDifference": error})
    predictions = defaultdict(list)
    for row in read("current_identifier_flagged_exclusion_predictions.csv"):
        predictions[row["Split"], row["Stage"], row["Threshold"]].append(row)
    flagged = {int(r["Index"]) for r in read(ROOT / "supplementary_audits/identity_sources/current_identifier_flagged.csv")}
    assert len(flagged) == 33
    for row in read("current_identifier_flagged_exclusion_cohort_counts.csv"):
        for part in ("train", "test"):
            full = official[row["Split"], part]
            retained = [i for i in full if i not in flagged]
            assert row[f"Original{part.title()}IndexSHA256"] == index_hash(full)
            assert row[f"Retained{part.title()}IndexSHA256"] == index_hash(retained)
            assert int(row[f"Retained{part.title()}N"]) == len(retained)
            membership_checks += 2
    estimates = {}
    for row in read("current_identifier_flagged_exclusion_auroc.csv"):
        key = row["Split"], row["Stage"], row["Threshold"]
        records = predictions[key]
        expected = official[row["Split"], "test"]
        if row["Stage"] != "original_all":
            expected = [i for i in expected if i not in flagged]
        assert [int(r["Index"]) for r in records] == expected
        assert all(int(r["Y"]) == original[int(row["Threshold"])][int(r["Index"])] for r in records)
        assert all(r["CAS"] == benchmark[int(r["Index"])]["CAS"] for r in records)
        membership_checks += 1
        value = auc([int(r["Y"]) for r in records], [float(r["Score"]) for r in records])
        error = abs(value - float(row["AUROC"]))
        assert error < 1e-12
        estimates[key] = value
        results.append({"Analysis": "identifier_flagged_exclusion", "Quantity": "AUROC", "MaxAbsDifference": error})
    for row in read("current_identifier_flagged_exclusion_contrasts.csv"):
        a, b, c = row["Split"], row["Stage"], row["Severe"]
        value = estimates[a, b, c] - estimates[a, b, "100"]
        error = abs(value - float(row["Gap"]))
        assert error < 1e-12
        results.append({"Analysis": "identifier_flagged_exclusion", "Quantity": "Gap", "MaxAbsDifference": error})
    report = {"Status": "passed", "Method": "independent explicit positive-negative pair enumeration with half-credit for ties",
              "Checks": len(results), "MaximumAbsoluteDifference": max(r["MaxAbsDifference"] for r in results),
              "MembershipHashAndPredictionOrderChecks": membership_checks,
              "Counts": {f"{a}/{q}": sum(r["Analysis"] == a and r["Quantity"] == q for r in results)
                         for a in ("uniform_D_exclude", "identifier_flagged_exclusion") for q in ("AUROC", "Gap")}}
    (OUT / "DERIVED_OUTPUT_VERIFICATION.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
