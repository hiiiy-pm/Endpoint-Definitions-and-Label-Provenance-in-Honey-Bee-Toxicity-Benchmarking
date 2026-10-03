"""Verify revision inputs, caches, model scores and label-audit outputs.

Run from any directory with the pinned analysis Python. Nothing is downloaded.
Cache negative tests and regeneration run exclusively in a temporary directory.
"""
from __future__ import annotations

import argparse
import ast
import gc
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import RobustScaler
from sklearn.svm import SVC

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
import _common as c
from _cache_validation import (CacheValidationError, audit_cache_from_raw,
                               runtime_environment, validate_cache)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def definitions(script: str, names: tuple[str, ...]) -> dict:
    """Load only named definitions to test actual rules without running analyses."""
    tree = ast.parse((ROOT / "code" / script).read_text(encoding="utf-8-sig"))
    selected = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    require(len(selected) == len(names), f"Missing test definitions in {script}")
    namespace = {"np": np, "pd": pd}
    exec(compile(ast.Module(body=selected, type_ignores=[]), script, "exec"), namespace)
    return namespace


def boundary_checks() -> dict:
    old = definitions("11_qualifier_propagation_audit.py", ("binary_label",))["binary_label"]
    new = definitions("12_label_rule_audit.py", ("group_label",))["group_label"]
    require(old(np.array([11., 12.])) == 0, "Historical overlapping =11 boundary changed")
    require(old(np.array([11.])) == 1, "Historical all-11 branch changed")
    cases = {
        "exact_11": ([11.], [11.], 1),
        "greater_than_11": ([np.nextafter(11., np.inf)], [np.inf], 0),
        "greater_or_equal_11": ([11.], [np.inf], -1),
        "exact_11_and_12": ([11., 12.], [11., 12.], -1),
        "greater_11_and_exact_12": ([np.nextafter(11., np.inf), 12.], [np.inf, 12.], 0),
    }
    for name, (lo, hi, expected) in cases.items():
        g = pd.DataFrame({"value": [11.] * len(lo), "lower": lo, "upper": hi})
        require(new(g, "C", 11.) == expected, f"Corrected boundary case failed: {name}")
    require(new(pd.DataFrame({"value": [11., 12.]}), "B", 11.) == -1,
            "Numeric unanimity must distinguish exact 11 from values >11")
    return {"historical_rule_cases": 2, "corrected_rule_cases": 6, "status": "passed"}


def negative_cache_checks(study) -> list[dict]:
    results = []
    with tempfile.TemporaryDirectory(prefix="apistox_cache_validation_") as temporary:
        base = Path(temporary)
        require(base.resolve().is_relative_to(Path(tempfile.gettempdir()).resolve()),
                "Temporary checks must stay inside the system temporary directory")

        def rejected(name, mutation, *, builder=False):
            directory = base / name
            shutil.copytree(c.CACHE, directory)
            mutation(directory)
            try:
                if builder:
                    c.build_or_load_cache(cache_dir=directory)
                else:
                    validate_cache(study, directory)
            except (CacheValidationError, ValueError) as error:
                results.append({"case": name, "status": "rejected", "reason": str(error)})
            else:
                raise AssertionError(f"Invalid cache was accepted: {name}")

        def metadata_change(directory, key, value):
            path = directory / "CACHE_METADATA.json"
            data = json.loads(path.read_text(encoding="utf-8"))
            data[key] = value
            path.write_text(json.dumps(data), encoding="utf-8")

        def swap_processed(directory):
            path = directory / "processed_data.csv"
            df = pd.read_csv(path)
            df.iloc[[0, 1]] = df.iloc[[1, 0]].to_numpy()
            df.to_csv(path, index=False)

        def bad_kernel(directory, kind):
            path = directory / "kernel_ecfp.npy"
            array = np.load(path)
            if kind == "nonfinite":
                array[0, 1] = np.nan
            elif kind == "asymmetry":
                array[0, 1] += .01
            elif kind == "permutation":
                idx = np.arange(len(array))
                idx[:2] = [1, 0]
                array = array[np.ix_(idx, idx)]
            elif kind == "shape":
                array = array[:-1, :-1]
            np.save(path, array)

        rejected("source_hash", lambda d: metadata_change(d, "source_sha256", "0" * 64))
        rejected("compound_row_order", swap_processed)
        rejected("nonfinite_kernel", lambda d: bad_kernel(d, "nonfinite"))
        rejected("asymmetric_kernel", lambda d: bad_kernel(d, "asymmetry"))
        rejected("reordered_symmetric_kernel", lambda d: bad_kernel(d, "permutation"))
        rejected("kernel_shape", lambda d: bad_kernel(d, "shape"))
        rejected("partial_cache", lambda d: (d / "kernel_maccs.npy").unlink(), builder=True)
        meta = json.loads((c.CACHE / "CACHE_METADATA.json").read_text(encoding="utf-8"))
        if meta.get("schema_version", 1) == 1:
            rejected("legacy_without_validation", lambda d: (d / "CACHE_VALIDATION.json").unlink())

        # Exercise genuine regeneration metadata, not a retrofitted legacy manifest.
        rebuilt = base / "new_cache"
        fresh = c.build_or_load_cache(force=True, cache_dir=rebuilt)
        metadata = json.loads((rebuilt / "CACHE_METADATA.json").read_text(encoding="utf-8"))
        require(metadata["schema_version"] == 2 and fresh["validation"]["status"] == "passed",
                "Fresh-cache metadata validation failed")
        del fresh  # Close mmap file handles before TemporaryDirectory cleanup on Windows.
        results.append({"case": "new_cache_v2_round_trip", "status": "passed"})
        for field in ("ordered_smiles_sha256", "ordered_identity_sha256", "input_sha256", "environment", "parameters"):
            damaged = dict(metadata)
            damaged[field] = "deliberately_invalid"
            (rebuilt / "CACHE_METADATA.json").write_text(json.dumps(damaged), encoding="utf-8")
            try:
                validate_cache(study, rebuilt)
            except CacheValidationError as error:
                results.append({"case": f"v2_{field}", "status": "rejected", "reason": str(error)})
            else:
                raise AssertionError(f"Invalid v2 metadata was accepted: {field}")

        split_dir = base / "duplicate_split"
        shutil.copytree(c.SPLIT_DIR, split_dir)
        path = split_dir / "random_train.csv"
        rows = pd.read_csv(path)
        rows.iloc[1] = rows.iloc[0]
        rows.to_csv(path, index=False)
        try:
            c.load_study_data(split_dir=split_dir)
        except AssertionError as error:
            results.append({"case": "duplicate_split_compound", "status": "rejected", "reason": str(error)})
        else:
            raise AssertionError("Duplicate split compound was accepted")
        gc.collect()
    return results


def fit_score(fresh, representation, train, test, labels):
    if representation == "Descriptors12":
        x = fresh["descriptors"]
        scaler = RobustScaler().fit(x[train])
        model = SVC(C=1, kernel="rbf", gamma="scale", class_weight="balanced")
        return model.fit(scaler.transform(x[train]), labels[train]).decision_function(scaler.transform(x[test]))
    kernel = fresh["kernels"][representation]
    model = SVC(C=1, kernel="precomputed", class_weight="balanced")
    return model.fit(kernel[np.ix_(train, train)], labels[train]).decision_function(kernel[np.ix_(test, train)])


def primary_score_checks(study, fresh, output) -> dict:
    saved = pd.read_csv(c.RESULTS / "primary/02_test_predictions.csv", float_precision="round_trip")
    rows = []
    for split, (train, test) in study.splits.items():
        for rep in c.REPRESENTATIONS:
            for threshold in c.THRESHOLDS:
                group = saved[(saved.Split == split) & (saved.Representation == rep)
                              & (saved.Threshold == threshold)].sort_values("TestOrder")
                require(np.array_equal(group.Index, test), "Primary prediction identity/order mismatch")
                require(np.array_equal(group.Y, study.y[threshold][test]), "Primary prediction labels mismatch")
                score = fit_score(fresh, rep, train, test, study.y[threshold])
                difference = float(np.max(np.abs(score - group.Score.to_numpy())))
                require(difference < 1e-8, f"Raw-SMILES refit differs: {split}/{rep}/{threshold}: {difference}")
                rows.append({"Split": split, "Representation": rep, "Threshold": threshold,
                             "n_test": len(test), "max_score_abs_difference": difference,
                             "AUROC": roc_auc_score(study.y[threshold][test], score)})
    pd.DataFrame(rows).to_csv(output / "primary_raw_smiles_refit.csv", index=False)
    return {"models": len(rows), "scores": sum(row["n_test"] for row in rows),
            "max_score_abs_difference": max(row["max_score_abs_difference"] for row in rows)}


def ecotox_checks(study, output) -> dict:
    # A second vectorized reproduction, separate from scripts 11/12 and their outputs.
    raw_path = ROOT / "output/external_data_audit/cache/upstream_ecotox.csv"
    raw = pd.read_csv(raw_path, sep="|", dtype=str)
    raw.columns = raw.columns.str.strip()
    units = {"AI ug/org": 1., "AI ug/org/d": 1., "ug/bee": 1., "ug/org": 1., "ug/org/d": 1.,
             "AI ng/org": .001, "AI ng/org/d": .001, "ng/org": .001,
             "AI mg/org": 1000., "mg/bee": 1000., "mg/org": 1000., "pg/org": .000001}
    routes = {"Diet, unspecified": "Oral", "Drinking water": "Oral", "Food": "Oral",
              "Dermal": "Contact", "Direct application": "Contact", "Topical, general": "Contact",
              "Multiple routes between application groups": "Other", "Oral via capsule": "Other",
              "Spray, unspecified": "Other", "Environmental, unspecified": "Other"}
    value = raw["Observed Response Mean"].str.strip()
    factor = raw["Observed Response Units"].str.strip().map(units)
    eligible = factor.notna() & value.ne("NR")
    records = pd.DataFrame({"CAS_numeric": raw.loc[eligible, "CAS Number"].astype(float).astype(int).astype(str),
                            "value": value[eligible].str.replace("/", "", regex=False).astype(float) * factor[eligible],
                            "route": raw.loc[eligible, "Exposure Type"].str.strip().map(routes)})
    records["CAS"] = records.CAS_numeric.map(lambda s: f"{s[:-3]}-{s[-3:-1]}-{s[-1]}")
    groups = records.groupby(["CAS", "route"]).value.agg(["min", "max", "median"]).reset_index()
    groups["rep_label"] = np.select([groups["max"] <= 11, groups["min"] >= 11], [1, 0], default=-1)
    groups["rep_level"] = np.select([groups["median"] > 100, groups["median"] > 1], [0, 1], default=2)
    groups["route_order"] = groups.route.map({"Contact": 0, "Oral": 1, "Other": 2})
    selected = groups[groups.rep_label != -1].sort_values(["CAS", "median", "route_order"]).drop_duplicates("CAS")
    joined = study.df[study.df.source.eq("ECOTOX")].merge(selected, on="CAS", how="left", validate="one_to_one")
    joined["match"] = (joined.label.eq(joined.rep_label) & joined.ppdb_level.eq(joined.rep_level)
                       & joined.toxicity_type.eq(joined.route))
    joined.to_csv(output / "ecotox_independent_reproduction.csv", index=False)
    require(len(joined) == 441 and joined.match.all(), "Independent ECOTOX reproduction is not 441/441")
    summary = json.loads((c.RESULTS / "tier_boundary/QUALIFIER_PROPAGATION_SUMMARY.json").read_text())
    require(summary["replication_label_and_level_match"] == 441, "Script 11 summary is not 441/441")
    return {"compounds": 441, "label_level_route_matches": 441, "raw_export_sha256": c.sha256_file(raw_path)}


def index_hash(indices) -> str:
    return hashlib.sha256(np.asarray(indices, dtype="<i8").tobytes()).hexdigest()


def crossover_labels(study) -> dict:
    """Reconstruct correction/eligibility independently from the exported rule ledger."""
    audit = pd.read_csv(c.RESULTS / "label_audit/02_compound_labels_by_rule.csv")
    labels = {}
    for scenario in ("primary", "conservative", "adult48"):
        labels[scenario] = {}
        for threshold in c.THRESHOLDS:
            rule = "C" if scenario == "conservative" or threshold == 11 else "D"
            filt = "adult48" if scenario == "adult48" else "all"
            group = audit[audit["filter"].eq(filt) & audit.routes.eq("eligible")
                          & audit.threshold.eq(threshold) & audit.rule.eq(rule)]
            require(len(group) == 441 and group.Index.is_unique, "Correction ledger is incomplete")
            y = study.y[threshold].astype(int).copy()
            y[group.Index.to_numpy(int)] = np.maximum(-1, group.label.to_numpy(int))
            labels[scenario][threshold] = y
    external = ROOT / "output/external_data_audit"
    concordance = pd.read_csv(external / "same_route_label_concordance.csv")
    concordance = concordance[concordance.source.isin(["OpenFoodTox", "PLOS2022"])]
    structures = pd.read_csv(external / "apistox_structures.csv")
    identity = dict(zip(structures.inchikey, structures.SMILES))
    row = {smiles: i for i, smiles in enumerate(study.df.SMILES)}
    labels["regulatory"], labels["combined"] = {}, {}
    for threshold in c.THRESHOLDS:
        y = study.y[threshold].astype(int).copy()
        combined = labels["primary"][threshold].copy()
        for key, group in concordance.groupby("inchikey"):
            i = row[identity[key]]
            outcomes = group[f"external_y_{threshold}"].astype(int).to_numpy()
            value = int(outcomes[0]) if (outcomes >= 0).all() and len(set(outcomes)) == 1 else -1
            y[i] = value
            combined[i] = (value if study.df.source.iloc[i] != "ECOTOX" or combined[i] == value else -1)
        labels["regulatory"][threshold] = y
        labels["combined"][threshold] = combined
    return labels


def crossover_checks(study, fresh) -> dict:
    directory = c.RESULTS / "label_audit"
    predictions = pd.read_csv(directory / "09_crossover_predictions.csv.gz", float_precision="round_trip")
    metrics = pd.read_csv(directory / "06_crossover_auroc.csv", float_precision="round_trip")
    gaps = pd.read_csv(directory / "07_crossover_gap.csv", float_precision="round_trip")
    keys = ["Scenario", "Unresolved", "Split", "Representation", "Threshold"]
    require(len(metrics) == 2700 and len(gaps) == 2700, "Crossover reporting grid is incomplete")
    expected_metrics = {}
    corrections = crossover_labels(study)
    max_difference = 0.
    max_score_difference = 0.
    exclude_refits = 0
    for key, group in predictions.groupby(keys, sort=False):
        scenario, how, split, rep, threshold = key
        train, official_test = study.splits[split]
        corrected = corrections[scenario]
        eligible = np.logical_and.reduce([corrected[t] >= 0 for t in c.THRESHOLDS])
        cohort = eligible if how == "exclude" else np.ones(len(study.df), dtype=bool)
        train = train[cohort[train]]
        test = group.Index.to_numpy(int)
        require(np.array_equal(test, official_test[cohort[official_test]]), "Crossover test identity/cohort mismatch")
        require(np.array_equal(group.OriginalLabel, study.y[threshold][test]), "Crossover original label mismatch")
        require(np.array_equal(group.CAS, study.df.CAS.to_numpy()[test]), "Crossover CAS identity mismatch")
        require(group.test_index_sha256.nunique() == 1 and group.test_index_sha256.iloc[0] == index_hash(test),
                "Crossover test membership hash mismatch")
        require(group.train_index_sha256.nunique() == 1 and group.train_index_sha256.iloc[0] == index_hash(train),
                "Crossover training membership hash mismatch")
        corrected_y = corrected[threshold].copy()
        if how != "exclude":
            replacement = 0 if how == "negative" and threshold == 100 else study.y[threshold]
            corrected_y = np.where(corrected_y < 0, replacement, corrected_y)
        require(np.array_equal(group.CorrectedLabel, corrected_y[test]), "Crossover corrected labels mismatch")
        for train_kind, score_column in (("original", "OriginalTrainScore"), ("corrected", "CorrectedTrainScore")):
            for eval_kind, label_column in (("original", "OriginalLabel"), ("corrected", "CorrectedLabel")):
                expected_metrics[(*key, train_kind, eval_kind)] = roc_auc_score(group[label_column], group[score_column])
        # Refit all exclude cells from independently reconstructed retained identities.
        if how == "exclude":
            for labels, column in ((study.y[threshold], "OriginalTrainScore"), (corrected_y, "CorrectedTrainScore")):
                score = fit_score(fresh, rep, train, test, labels)
                difference = float(np.max(np.abs(score - group[column].to_numpy())))
                max_score_difference = max(max_score_difference, difference)
                require(difference < 1e-8, f"Retained-cohort crossover refit differs for {key}/{column}")
            exclude_refits += 1
    require(len(expected_metrics) == len(metrics), "Duplicate or missing crossover model cells")
    for row in metrics.itertuples(index=False):
        key = tuple(getattr(row, key) for key in keys)
        actual = expected_metrics[(*key, row.TrainLabels, row.EvalLabels)]
        max_difference = max(max_difference, abs(actual - row.AUROC))
        require(abs(actual - row.AUROC) < 1e-12, "Crossover AUROC does not match saved predictions")
        require(row.B_valid + row.B_degenerate == row.B_requested, "Invalid AUROC bootstrap accounting")
    for row in gaps.itertuples(index=False):
        key = (row.Scenario, row.Unresolved, row.Split, row.Representation)
        def gap(a, b):
            return (expected_metrics[(*key, row.Severe, a, b)]
                    - expected_metrics[(*key, 100, a, b)])
        if row.TrainLabels.startswith("change:"):
            train_kind = row.TrainLabels.split(":")[1].split("-vs-")[0]
            expected = gap(train_kind, "corrected") - gap("original", "original")
        else:
            expected = gap(row.TrainLabels, row.EvalLabels)
        require(abs(expected - row.Gap) < 1e-12, "Crossover gap is inconsistent with AUROCs")
        require(row.B_valid + row.B_degenerate == row.B_requested, "Invalid gap bootstrap accounting")
    return {"prediction_rows": len(predictions), "AUROC_cells": len(metrics), "gap_cells": len(gaps),
            "max_AUROC_abs_difference": max_difference, "exclude_model_cells": exclude_refits,
            "exclude_models_refitted": 2 * exclude_refits,
            "exclude_refit_max_score_abs_difference": max_score_difference,
            "scope": "Correction ledgers, cohort identities, all metric/gap algebra, bootstrap accounting, and both training-label models in all exclude cells"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-only", action="store_true", help="Only recompute/validate the existing cache")
    parser.add_argument("--output", type=Path, default=ROOT / "results/reproducibility")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    study = c.load_study_data()
    report, fresh = audit_cache_from_raw(study, c.CACHE)
    result = {"status": "passed", "environment": runtime_environment(), "cache": report}
    print("Raw-SMILES cache comparison passed", flush=True)
    if not args.cache_only:
        result["negative_paths"] = negative_cache_checks(study)
        result["boundary_rules"] = boundary_checks()
        result["primary_refits"] = primary_score_checks(study, fresh, args.output)
        print("Cache guards and 45 primary-model refits passed", flush=True)
        result["ecotox_replication"] = ecotox_checks(study, args.output)
        result["crossover"] = crossover_checks(study, fresh)
    (args.output / "revision_verification.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in result.items() if key != "cache"}, indent=2))


if __name__ == "__main__":
    main()
