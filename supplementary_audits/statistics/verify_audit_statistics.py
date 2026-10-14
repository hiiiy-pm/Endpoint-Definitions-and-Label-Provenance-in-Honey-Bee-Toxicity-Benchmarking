"""Independent, standard-library reconstruction of label-rule sensitivity.

Frozen inputs are read without modification. Outputs are confined to this
directory. The Benomyl operator change is a
counterfactual replay, not an edit of the archived ECOTOX export.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import median

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
THRESHOLDS = (100, 11, 1)
EPS = 1e-9
UNIT_FACTOR = {"AI ug/org": 1, "AI ug/org/d": 1, "ug/bee": 1,
               "ug/org": 1, "ug/org/d": 1, "AI ng/org": 1e-3,
               "AI ng/org/d": 1e-3, "ng/org": 1e-3, "AI mg/org": 1e3,
               "mg/bee": 1e3, "mg/org": 1e3, "pg/org": 1e-6}
ROUTE = {"Diet, unspecified": "Oral", "Drinking water": "Oral", "Food": "Oral",
         "Dermal": "Contact", "Direct application": "Contact", "Topical, general": "Contact",
         "Multiple routes between application groups": "Other", "Oral via capsule": "Other",
         "Spray, unspecified": "Other", "Environmental, unspecified": "Other"}


def read(relative, delimiter=","):
    with (ROOT / relative).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_csv(name, rows):
    assert rows, name
    with (OUT / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tier(row):
    level = int(row["ppdb_level"])
    return 0 if level == 0 else (3 if level == 2 else 1 + int(row["label"]))


def binary(values):
    return 1 if max(values) <= 11 else (0 if min(values) >= 11 else -1)


def record_interval(row):
    op, value = row["operator"], row["value"]
    lo = value if op in ("", "~", ">=") else (value * (1 + EPS) if op == ">" else 0.0)
    hi = value if op in ("", "~", "<=") else (value * (1 - EPS) if op == "<" else math.inf)
    return lo, hi


def group_label(rows, rule, threshold):
    if not rows:
        return -2
    values = [r["value"] for r in rows]
    if rule == "A":
        return int(median(values) <= threshold)
    if rule == "B":
        return 1 if max(values) <= threshold else (0 if min(values) > threshold else -1)
    lower, upper = zip(*(record_interval(r) for r in rows))
    if rule == "C":
        return 1 if max(upper) <= threshold else (0 if min(lower) > threshold else -1)
    assert rule == "D"
    return 1 if median(upper) <= threshold else (0 if median(lower) > threshold else -1)


def route_min(labels):
    labels = [y for y in labels if y != -2]
    return -2 if not labels else (1 if 1 in labels else (0 if set(labels) == {0} else -1))


def main():
    benchmark = read("data/raw/apistox.csv")
    assert len(benchmark) == 1035
    original = {t: [int(tier(r) >= {100: 1, 11: 2, 1: 3}[t]) for r in benchmark]
                for t in THRESHOLDS}
    audit = read("results/label_audit/02_compound_labels_by_rule.csv")
    lookup = {(int(r["Index"]), r["filter"], r["routes"], int(r["threshold"]), r["rule"]): int(r["label"])
              for r in audit}
    assert len(lookup) == len(audit)
    for row in audit:
        i, t = int(row["Index"]), int(row["threshold"])
        assert row["CAS"] == benchmark[i]["CAS"] and int(row["original"]) == original[t][i]

    def scenario(all_d):
        result = {t: original[t].copy() for t in THRESHOLDS}
        for t in THRESHOLDS:
            rule = "D" if all_d or t != 11 else "C"
            selected = [r for r in audit if r["filter"] == "all" and r["routes"] == "eligible"
                        and int(r["threshold"]) == t and r["rule"] == rule]
            assert len(selected) == 441
            for row in selected:
                result[t][int(row["Index"])] = max(-1, int(row["label"]))
        return result

    mixed, uniform = scenario(False), scenario(True)
    raw_differences, vector_rows = [], []
    for t in THRESHOLDS:
        for i in range(len(benchmark)):
            if mixed[t][i] != uniform[t][i]:
                raw_differences.append({"Index": i, "CAS": benchmark[i]["CAS"], "Name": benchmark[i]["name"],
                                       "Tier": tier(benchmark[i]), "Threshold": t, "Original": original[t][i],
                                       "Mixed_DCD": mixed[t][i], "Uniform_DDD": uniform[t][i]})
        a = [o if y == -1 else y for y, o in zip(mixed[t], original[t])]
        b = [o if y == -1 else y for y, o in zip(uniform[t], original[t])]
        vector_rows.append({"Threshold": t, "Length": len(a), "KeepDifferences": sum(x != y for x, y in zip(a, b)),
                            "MixedKeepPositive": sum(a), "UniformKeepPositive": sum(b),
                            "ChangedFromOriginal": sum(x != y for x, y in zip(a, original[t])),
                            "MixedKeepUint8SHA256": hashlib.sha256(bytes(a)).hexdigest(),
                            "UniformKeepUint8SHA256": hashlib.sha256(bytes(b)).hexdigest()})
        assert a == b
    assert {(r["Index"], r["Threshold"]) for r in raw_differences} == {(861, 11), (953, 11)}
    km = [all(mixed[t][i] != -1 for t in THRESHOLDS) for i in range(len(benchmark))]
    ku = [all(uniform[t][i] != -1 for t in THRESHOLDS) for i in range(len(benchmark))]
    additional = [i for i in range(len(benchmark)) if ku[i] and not km[i]]
    assert (sum(km), sum(ku), additional) == (913, 915, [861, 953])
    assert not any(a and not b for a, b in zip(km, ku))
    assert all(mixed[t][i] == uniform[t][i] for i in range(len(benchmark)) if km[i] for t in THRESHOLDS)
    split_rows, additions = [], []
    split_members = {}
    for split in ("random", "maxmin", "time"):
        for partition in ("train", "test"):
            smiles = {r["SMILES"] for r in read(f"data/official_splits/{split}_{partition}.csv")}
            indices = [i for i, r in enumerate(benchmark) if r["SMILES"] in smiles]
            split_members[split, partition] = indices
            extra = [i for i in additional if i in indices]
            split_rows.append({"Split": split, "Partition": partition, "MixedExcludeN": sum(km[i] for i in indices),
                               "UniformExcludeN": sum(ku[i] for i in indices), "AddedN": len(extra),
                               "AddedIndices": ";".join(map(str, extra))})
            for i in extra:
                additions.append({"Split": split, "Partition": partition, "Index": i,
                                  "CAS": benchmark[i]["CAS"], "Name": benchmark[i]["name"],
                                  **{f"Y{t}": uniform[t][i] for t in THRESHOLDS}})
    write_csv("all_d_vector_comparison.csv", vector_rows)
    write_csv("all_d_raw_label_differences.csv", raw_differences)
    write_csv("all_d_exclusion_membership.csv", split_rows)
    write_csv("all_d_exclusion_additions.csv", additions)

    # Reconstruct the logging/filtering distinction without importing upstream code.
    raw = [{k.strip(): (v or "").strip() for k, v in row.items() if k is not None}
           for row in read("output/external_data_audit/cache/upstream_ecotox.csv", "|")]
    records, groups = [], defaultdict(list)
    for ordinal, row in enumerate(raw):
        unit = row["Observed Response Units"]
        if unit not in UNIT_FACTOR or row["Observed Response Mean"] == "NR":
            continue
        digits = str(int(float(row["CAS Number"])))
        record = {"SourceRow0": ordinal, "CAS": f"{digits[:-3]}-{digits[-3:-1]}-{digits[-1]}",
                  "route": ROUTE.get(row["Exposure Type"]), "operator": row["Observed Response Mean Op"],
                  "value": float(row["Observed Response Mean"].replace("/", "")) * UNIT_FACTOR[unit],
                  "unit": unit, "stage": row["Organism Lifestage"],
                  "duration": row["Observed Duration (Days)"], "duration_op": row["Observed Duration Op (Days)"],
                  "duration_unit": row["Observed Duration Units (Days)"], "conc_type": row["Conc 1 Type (Author)"],
                  "reference": row["Reference Number"]}
        record["adult48"] = (record["stage"] == "Adult" and record["duration"] == "2"
                             and record["duration_op"] == "" and record["duration_unit"] == "Day(s)"
                             and not unit.endswith("/d"))
        record["adult48ai"] = record["adult48"] and record["conc_type"] == "Active ingredient"
        records.append(record)
        if record["route"]:
            groups[record["CAS"], record["route"]].append(record)
    removed_groups = {key for key, values in groups.items() if binary([r["value"] for r in values]) == -1}
    logged_cas = {key[0] for key in removed_groups}
    logged = [r for r in records if r["CAS"] in logged_cas]
    removed = [r for r in logged if (r["CAS"], r["route"]) in removed_groups]
    retained = [r for r in logged if (r["CAS"], r["route"]) not in removed_groups]
    eco = {r["CAS"]: r for r in benchmark if r["source"] == "ECOTOX"}
    summary_log = {"RawRecords": len(raw), "UnitValueEligibleRecords": len(records),
                   "CASRouteGroups": len(groups), "RemovedCASRouteGroups": len(removed_groups),
                   "LoggedCASIdentities": len(logged_cas), "LoggedRecords": len(logged),
                   "ActualRemovedGroupRecords": len(removed), "LoggedRetainedGroupRecords": len(retained),
                   "LoggedRetainedGroupCAS": len({r["CAS"] for r in retained}),
                   "LoggedCASInFinalBenchmark": len(logged_cas & {r["CAS"] for r in benchmark}),
                   "LoggedRetainedCASInFinalECOTOX": len({r["CAS"] for r in retained} & set(eco)),
                   "LoggedRetainedFinalDeterminingRecords": sum(r["CAS"] in eco and r["route"] == eco[r["CAS"]]["toxicity_type"] for r in retained)}
    assert (len(logged), len(removed), len(retained)) == (745, 536, 209)
    assert summary_log["LoggedRetainedFinalDeterminingRecords"] == 119
    write_csv("exclusion_log_summary.csv", [summary_log])
    write_csv("exclusion_log_record_replay.csv", [{"SourceRow0": r["SourceRow0"], "CAS": r["CAS"], "Route": r["route"],
               "ActualGroupExcluded": (r["CAS"], r["route"]) in removed_groups,
               "FinalECOTOXDeterminingGroup": r["CAS"] in eco and r["route"] == eco[r["CAS"]]["toxicity_type"]} for r in logged])

    # Counterfactual Benomyl replay: restore the reported >10 sign while retaining
    # the archived numerical values and historical route eligibility.
    bi = next(i for i, r in enumerate(benchmark) if r["CAS"] == "17804-35-2")
    assert bi == 111
    benomyl = [r for r in records if r["CAS"] == benchmark[bi]["CAS"]]
    disputed = [r for r in benomyl if r["route"] == "Oral" and r["operator"] == "<"
                and r["value"] == 10 and r["reference"] == "35461"]
    assert len(disputed) == 1 and disputed[0]["SourceRow0"] == 1251
    changed_records = [dict(r, operator=">") if r["SourceRow0"] == 1251 else dict(r) for r in benomyl]
    eligible_routes = {route for cas, route in groups if cas == benchmark[bi]["CAS"] and (cas, route) not in removed_groups}

    def ben_label(rows, filt, routes, rule, t):
        selected = rows if filt == "all" else [r for r in rows if r[filt]]
        route_set = ({benchmark[bi]["toxicity_type"]} if routes == "fixed" else
                     eligible_routes if routes == "eligible" else {r["route"] for r in rows})
        labels = [group_label([r for r in selected if r["route"] == route], rule, t) for route in sorted(route_set)]
        return labels[0] if routes == "fixed" else route_min(labels)

    replay = []
    for filt in ("all", "adult48", "adult48ai"):
        for routes in ("fixed", "eligible", "all"):
            for t in THRESHOLDS:
                for rule in ("A", "B", "C", "D"):
                    before = ben_label(benomyl, filt, routes, rule, t)
                    after = ben_label(changed_records, filt, routes, rule, t)
                    assert before == lookup[bi, filt, routes, t, rule]
                    replay.append({"Index": bi, "CAS": benchmark[bi]["CAS"], "Filter": filt, "Routes": routes,
                                   "Threshold": t, "Rule": rule, "ArchivedLabel": before,
                                   "ReplayedGreaterThan10Label": after, "Changed": before != after})
    write_csv("benomyl_rule_replay.csv", replay)
    structure = read("output/external_data_audit/apistox_structures.csv")
    ben_key = next(r["inchikey"] for r in structure if r["SMILES"] == benchmark[bi]["SMILES"])
    regulatory = [r for r in read("output/external_data_audit/same_route_label_concordance.csv")
                  if r["inchikey"] == ben_key and r["source"] in ("OpenFoodTox", "PLOS2022")]

    def ben_scenario(rows, name):
        filt = "adult48" if name == "adult48" else "all"
        base = {t: max(-1, ben_label(rows, filt, "eligible", "C" if name == "conservative" or t == 11 else "D", t))
                for t in THRESHOLDS}
        if name in ("regulatory", "combined"):
            for t in THRESHOLDS:
                if not regulatory:
                    base[t] = original[t][bi] if name == "regulatory" else base[t]
                else:
                    vals = {int(r[f"external_y_{t}"]) for r in regulatory}
                    reg = next(iter(vals)) if len(vals) == 1 and -1 not in vals else -1
                    base[t] = reg if name == "regulatory" or reg == base[t] else -1
        return base

    def handle(labels, policy):
        retained = policy != "exclude" or all(y != -1 for y in labels.values())
        final = {t: (0 if policy == "negative" and t == 100 else original[t][bi]) if y == -1 and policy != "exclude" else y
                 for t, y in labels.items()}
        return retained, final

    ben_scenario_rows = []
    for name in ("primary", "conservative", "adult48", "regulatory", "combined"):
        before, after = ben_scenario(benomyl, name), ben_scenario(changed_records, name)
        for policy in ("keep", "exclude", "negative"):
            kb, lb = handle(before, policy)
            ka, la = handle(after, policy)
            ben_scenario_rows.append({"Scenario": name, "Policy": policy, "Index": bi, "CAS": benchmark[bi]["CAS"],
                "ArchivedRetained": kb, "ReplayedRetained": ka,
                **{f"ArchivedY{t}": lb[t] for t in THRESHOLDS}, **{f"ReplayedY{t}": la[t] for t in THRESHOLDS},
                "EffectiveInputChanged": kb != ka or (kb and ka and lb != la),
                "RandomMembership": "test", "MaxMinMembership": "train", "TimeMembership": "train"})
    assert all(not r["EffectiveInputChanged"] for r in ben_scenario_rows if r["Policy"] in ("keep", "exclude"))
    write_csv("benomyl_scenario_input_impact.csv", ben_scenario_rows)

    cid_groups = defaultdict(list)
    for i, row in enumerate(benchmark):
        cid_groups[row["CID"]].append((i, row))
    duplicates = {cid: rows for cid, rows in cid_groups.items() if len(rows) > 1}
    assert len(duplicates) == 5 and sum(map(len, duplicates.values())) == 10
    duplicate_rows = []
    for cid, rows in duplicates.items():
        for i, row in rows:
            duplicate_rows.append({"CID": cid, "Index": i, "CSVLine1": i + 2, "CAS": row["CAS"],
                                   "Name": row["name"], "SMILES": row["SMILES"],
                                   "Category": "cation_and_chloride_salt" if cid == "13836" else
                                               "metiram_zineb_representation_conflict" if cid == "5284484" else "different_structures",
                                   "Adjudication": "duplicate identifier confirmed; correct CID requires external verification"})
    write_csv("duplicate_cid_inventory.csv", duplicate_rows)

    inputs = ["data/raw/apistox.csv", "results/label_audit/02_compound_labels_by_rule.csv",
              "output/external_data_audit/cache/upstream_ecotox.csv", "output/external_data_audit/apistox_structures.csv",
              "output/external_data_audit/same_route_label_concordance.csv", "code/13_corrected_label_crossover.py"]
    inputs += [f"data/official_splits/{s}_{p}.csv" for s in ("random", "maxmin", "time") for p in ("train", "test")]
    result = {"Status": "passed", "Method": "independent standard-library reconstruction; no model fitting",
              "AllDKeepVectorsIdentical": True, "MixedExcludeN": sum(km), "UniformExcludeN": sum(ku),
              "AddedIndices": additional, "ExclusionLog": summary_log,
              "BenomylArchivedOperator": "<", "BenomylReplayedOperator": ">", "BenomylSourceRow0": 1251,
              "BenomylKeepExcludeInputsUnchanged": True,
              "DuplicateCIDGroups": len(duplicates), "DuplicateCIDRows": len(duplicate_rows),
              "InputSHA256": {p: digest(ROOT / p) for p in inputs},
              "Limits": ["Benomyl replay does not adjudicate the primary paper or edit the archive.",
                         "Identical keep inputs establish deterministic equivalence under the same fitting protocol.",
                         "All-D exclusion changes cohort membership; performance is evaluated separately by run_uniform_d_exclude.py."]}
    (OUT / "STATISTICS_VERIFICATION.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "InputSHA256"}, indent=2))


if __name__ == "__main__":
    main()
