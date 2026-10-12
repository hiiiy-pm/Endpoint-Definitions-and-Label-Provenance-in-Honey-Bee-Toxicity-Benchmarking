"""ECFP sensitivity to uniform interval-median aggregation with exclusion.

All outputs are separate from the frozen primary results. The original full-
cohort fits and the published D/C/D-exclude fits must reproduce saved scores
before any new result is exported. New scenario comparisons use common
scaffold draws on the union test cohort and are conditional on fitted models.
"""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import roc_auc_score
from sklearn.svm import SVC

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT / "code"))
from _common import load_study_data, multinomial_weights, percentile_ci, tie_aware_weighted_auc

B = 2000
SEED = 20260829
THRESHOLDS = (100, 11, 1)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    study = load_study_data()
    kernel = np.load(ROOT / "results/cache/kernel_ecfp.npy")
    scaffolds = pd.read_csv(ROOT / "results/cache/scaffolds.csv").sort_values("Index").Scaffold.to_numpy(str)
    original = study.y
    audit = pd.read_csv(ROOT / "results/label_audit/02_compound_labels_by_rule.csv")
    saved_primary = pd.read_csv(ROOT / "results/primary/02_test_predictions.csv")
    saved_crossover = pd.read_csv(ROOT / "results/label_audit/09_crossover_predictions.csv.gz")
    saved_crossover = saved_crossover[(saved_crossover.Scenario == "primary") &
                                      (saved_crossover.Unresolved == "exclude") &
                                      (saved_crossover.Representation == "ECFP")]
    labels, kept = {}, {}
    for scenario in ("mixed_DCD", "uniform_DDD"):
        labels[scenario] = {t: original[t].copy() for t in THRESHOLDS}
        for t in THRESHOLDS:
            rule = "C" if scenario == "mixed_DCD" and t == 11 else "D"
            rows = audit[(audit["filter"] == "all") & (audit.routes == "eligible") &
                         (audit.threshold == t) & (audit["rule"] == rule)]
            assert len(rows) == 441
            idx = rows.Index.to_numpy(int)
            assert np.array_equal(rows.original.to_numpy(int), original[t][idx])
            labels[scenario][t][idx] = np.maximum(rows.label.to_numpy(int), -1)
        kept[scenario] = np.logical_and.reduce([labels[scenario][t] != -1 for t in THRESHOLDS])
    assert kept["mixed_DCD"].sum() == 913 and kept["uniform_DDD"].sum() == 915
    assert set(np.flatnonzero(kept["uniform_DDD"] & ~kept["mixed_DCD"])) == {861, 953}

    def fit(tr, y):
        return SVC(C=1, kernel="precomputed", class_weight="balanced").fit(kernel[np.ix_(tr, tr)], y[tr])

    def predict(model, tr, te):
        return model.decision_function(kernel[np.ix_(te, tr)])

    checks = []
    for split, (tr, te) in study.splits.items():
        for t in THRESHOLDS:
            fitted = fit(tr, original[t])
            scores = predict(fitted, tr, te)
            saved = saved_primary[(saved_primary.Split == split) & (saved_primary.Representation == "ECFP") &
                                  (saved_primary.Threshold == t)].sort_values("TestOrder")
            assert np.array_equal(saved.Index.to_numpy(int), te)
            error = float(np.max(np.abs(scores - saved.Score.to_numpy(float))))
            assert error < 1e-9, (split, t, error)
            checks.append({"Check": "full_original", "Split": split, "Threshold": t, "TrainLabels": "original",
                           "NTest": len(te), "MaxAbsoluteScoreDifference": error})

    auc_rows, gap_rows, sensitivity_rows, prediction_rows, fit_rows = [], [], [], [], []
    for split_index, (split, (tr, te)) in enumerate(study.splits.items()):
        # The uniform-D test cohort is the union; mixed-D/C/D membership is a subset.
        union_te = te[kept["uniform_DDD"][te]]
        clusters, inverse = np.unique(scaffolds[union_te], return_inverse=True)
        rng = np.random.default_rng(SEED + 13007 + 97 * split_index)
        union_weights = multinomial_weights(len(clusters), B, rng)[:, inverse]
        scenario_gaps, scenario_gap_draws = {}, {}
        fit_cache = {}
        for scenario in ("mixed_DCD", "uniform_DDD"):
            trk, tek = tr[kept[scenario][tr]], te[kept[scenario][te]]
            union_mask = np.isin(union_te, tek)
            assert np.array_equal(union_te[union_mask], tek)
            weights = union_weights[:, union_mask]
            est, draws = {}, {}
            metadata = {"Scenario": scenario, "Split": split, "n_train": len(trk), "n_test": len(tek),
                        "train_index_sha256": hashlib.sha256(trk.astype("<i8").tobytes()).hexdigest(),
                        "test_index_sha256": hashlib.sha256(tek.astype("<i8").tobytes()).hexdigest()}
            for t in THRESHOLDS:
                scores = {}
                for train_labels in ("original", "provisional"):
                    y = original[t] if train_labels == "original" else labels[scenario][t]
                    key = (t, train_labels, tuple(trk), tuple(y[trk]))
                    reused = key in fit_cache
                    model = fit_cache[key] if reused else fit(trk, y)
                    fit_cache[key] = model
                    scores[train_labels] = predict(model, trk, tek)
                    fit_rows.append({**metadata, "Threshold": t, "TrainLabels": train_labels,
                                     "ReusedIdenticalTrainingFit": reused})
                    if scenario == "mixed_DCD":
                        saved = saved_crossover[(saved_crossover.Split == split) & (saved_crossover.Threshold == t)]
                        saved = saved.set_index("Index").loc[tek]
                        assert (saved.n_train == len(trk)).all() and (saved.n_test == len(tek)).all()
                        column = "OriginalTrainScore" if train_labels == "original" else "CorrectedTrainScore"
                        error = float(np.max(np.abs(scores[train_labels] - saved[column].to_numpy(float))))
                        assert error < 1e-9, (scenario, split, t, train_labels, error)
                        checks.append({"Check": "mixed_exclude_saved", "Split": split, "Threshold": t,
                                       "TrainLabels": train_labels, "NTest": len(tek),
                                       "MaxAbsoluteScoreDifference": error})
                for j, i in enumerate(tek):
                    prediction_rows.append({**metadata, "Threshold": t, "Index": int(i), "CAS": study.df.iloc[i].CAS,
                                            "OriginalLabel": int(original[t][i]),
                                            "ProvisionalLabel": int(labels[scenario][t][i]),
                                            "OriginalTrainScore": float(scores["original"][j]),
                                            "ProvisionalTrainScore": float(scores["provisional"][j])})
                for train_labels in ("original", "provisional"):
                    for eval_labels in ("original", "provisional"):
                        y = (original[t] if eval_labels == "original" else labels[scenario][t])[tek]
                        key = (t, train_labels, eval_labels)
                        est[key] = float(roc_auc_score(y, scores[train_labels]))
                        draws[key] = tie_aware_weighted_auc(y, scores[train_labels], weights)
                        lo, hi = percentile_ci(draws[key])
                        auc_rows.append({**metadata, "Threshold": t, "TrainLabels": train_labels,
                                         "EvalLabels": eval_labels, "n_pos": int(y.sum()), "AUROC": est[key],
                                         "CI_low": lo, "CI_high": hi, "B_requested": B,
                                         "B_valid": int(np.isfinite(draws[key]).sum())})
            for severe in (11, 1):
                for train_labels in ("original", "provisional"):
                    for eval_labels in ("original", "provisional"):
                        key = (severe, train_labels, eval_labels)
                        gap = est[key] - est[100, train_labels, eval_labels]
                        gd = draws[key] - draws[100, train_labels, eval_labels]
                        scenario_gaps[scenario, *key] = gap
                        scenario_gap_draws[scenario, *key] = gd
                        lo, hi = percentile_ci(gd)
                        gap_rows.append({**metadata, "Severe": severe, "TrainLabels": train_labels,
                                         "EvalLabels": eval_labels, "Gap": gap, "CI_low": lo, "CI_high": hi,
                                         "B_requested": B, "B_valid": int(np.isfinite(gd).sum())})
                key_p = (scenario, severe, "provisional", "provisional")
                key_o = (scenario, severe, "original", "original")
                effect = scenario_gaps[key_p] - scenario_gaps[key_o]
                effect_draw = scenario_gap_draws[key_p] - scenario_gap_draws[key_o]
                lo, hi = percentile_ci(effect_draw)
                gap_rows.append({**metadata, "Severe": severe, "TrainLabels": "change:provisional-vs-original",
                                 "EvalLabels": "change:provisional-vs-original", "Gap": effect,
                                 "CI_low": lo, "CI_high": hi, "B_requested": B,
                                 "B_valid": int(np.isfinite(effect_draw).sum())})
                scenario_gaps[scenario, severe, "change", "change"] = effect
                scenario_gap_draws[scenario, severe, "change", "change"] = effect_draw
        for severe in (11, 1):
            for arm in ("original", "provisional", "change"):
                a = ("mixed_DCD", severe, arm, arm)
                b = ("uniform_DDD", severe, arm, arm)
                difference = scenario_gaps[b] - scenario_gaps[a]
                difference_draws = scenario_gap_draws[b] - scenario_gap_draws[a]
                lo, hi = percentile_ci(difference_draws)
                sensitivity_rows.append({"Split": split, "Severe": severe, "Quantity": arm,
                                         "MixedEstimate": scenario_gaps[a], "UniformEstimate": scenario_gaps[b],
                                         "UniformMinusMixed": difference, "CI_low": lo, "CI_high": hi,
                                         "B_requested": B, "B_valid": int(np.isfinite(difference_draws).sum()),
                                         "UnionTestN": len(union_te), "UnionScaffolds": len(clusters)})
    assert len(checks) == 27
    for filename, rows in [("uniform_d_exclude_baseline_checks.csv", checks), ("uniform_d_exclude_auroc.csv", auc_rows),
                           ("uniform_d_exclude_gap.csv", gap_rows), ("uniform_d_exclude_scenario_differences.csv", sensitivity_rows),
                           ("uniform_d_exclude_predictions.csv", prediction_rows), ("uniform_d_exclude_fit_log.csv", fit_rows)]:
        pd.DataFrame(rows).to_csv(OUT / filename, index=False)
    meta = {"Status": "completed", "Representation": "ECFP", "BootstrapDraws": B,
            "BaselineChecks": len(checks), "MaxBaselineScoreDifference": max(r["MaxAbsoluteScoreDifference"] for r in checks),
            "TimeUniformFitsReused": sum(r["ReusedIdenticalTrainingFit"] for r in fit_rows if r["Split"] == "Time" and r["Scenario"] == "uniform_DDD"),
            "Resampling": "paired scaffold-cluster bootstrap on the union (uniform-D) test cohort; mixed cohort filtered within each draw",
            "InferenceScope": "conditional on fitted models and selected cohorts; exploratory, no multiplicity correction",
            "Interpretation": "scenario changes include membership changes; keep remains exactly input-equivalent",
            "Python": platform.python_version(), "Numpy": np.__version__, "Pandas": pd.__version__, "ScikitLearn": sklearn.__version__,
            "InputSHA256": {p: sha(ROOT / p) for p in ["results/cache/kernel_ecfp.npy", "results/cache/scaffolds.csv",
                 "results/label_audit/02_compound_labels_by_rule.csv", "results/primary/02_test_predictions.csv",
                 "results/label_audit/09_crossover_predictions.csv.gz"]}}
    (OUT / "UNIFORM_D_EXCLUDE_METADATA.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    q = pd.DataFrame(gap_rows)
    q = q[(q.Scenario == "uniform_DDD") & q.TrainLabels.str.startswith("change:")]
    print(q[["Split", "Severe", "Gap", "CI_low", "CI_high", "B_valid"]].to_string(index=False))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
