"""Predefined current-identifier flag exclusion sensitivity for fixed ECFP models.

The 33 flags are current identifier-consistency findings, not adjudicated
historical structure errors. Frozen structures, labels and official partitions
are preserved. Report all three splits, three endpoints, and two contrasts.
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
IDENTITY = ROOT / "supplementary_audits/identity_sources"
sys.path.insert(0, str(ROOT / "code"))
from _common import load_study_data, multinomial_weights, percentile_ci, simultaneous_centered_ci, tie_aware_weighted_auc

B = 5000
SEED = 20260829
THRESHOLDS = (100, 11, 1)
STAGES = ("original_all", "frozen_reduced", "refit_reduced")
FLAG_HASH = "344510081bd181496db1062b7e16c44a93746c6a3b15bbabb409636b30769e12"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    flag_path = IDENTITY / "current_identifier_flagged.csv"
    freeze = json.loads((IDENTITY / "CURRENT_IDENTIFIER_FLAGGED_FREEZE.json").read_text(encoding="utf-8"))
    assert sha(flag_path) == FLAG_HASH == freeze["sha256"]
    flagged = pd.read_csv(flag_path)
    assert len(flagged) == flagged.Index.nunique() == 33
    assert set(flagged.Index) == set(freeze["indices"])
    study = load_study_data()
    for row in flagged.itertuples(index=False):
        assert study.df.iloc[row.Index].CAS == row.CAS
    retain = np.ones(len(study.df), dtype=bool)
    retain[flagged.Index.to_numpy(int)] = False
    kernel = np.load(ROOT / "results/cache/kernel_ecfp.npy")
    scaffolds = pd.read_csv(ROOT / "results/cache/scaffolds.csv").sort_values("Index").Scaffold.to_numpy(str)
    saved = pd.read_csv(ROOT / "results/primary/02_test_predictions.csv")
    checks, counts, predictions, auc_rows, contrast_rows, increments = [], [], [], [], [], []
    for split_index, (split, (tr, te)) in enumerate(study.splits.items()):
        trk, tek = tr[retain[tr]], te[retain[te]]
        mask = retain[te]
        counts.append({"Split": split, "OriginalTrainN": len(tr), "RetainedTrainN": len(trk),
                       "ExcludedTrainN": len(tr) - len(trk), "OriginalTestN": len(te),
                       "RetainedTestN": len(tek), "ExcludedTestN": len(te) - len(tek),
                       "ExcludedTrainIndices": ";".join(map(str, tr[~retain[tr]])),
                       "ExcludedTestIndices": ";".join(map(str, te[~retain[te]])),
                       "OriginalTrainIndexSHA256": hashlib.sha256(tr.astype("<i8").tobytes()).hexdigest(),
                       "RetainedTrainIndexSHA256": hashlib.sha256(trk.astype("<i8").tobytes()).hexdigest(),
                       "OriginalTestIndexSHA256": hashlib.sha256(te.astype("<i8").tobytes()).hexdigest(),
                       "RetainedTestIndexSHA256": hashlib.sha256(tek.astype("<i8").tobytes()).hexdigest()})
        unique, inverse = np.unique(scaffolds[te], return_inverse=True)
        rng = np.random.default_rng(SEED + 18000 + 101 * split_index)
        full_weights = multinomial_weights(len(unique), B, rng)[:, inverse]
        estimates, draws = {}, {}
        for threshold in THRESHOLDS:
            y = study.y[threshold]
            source = saved[(saved.Split == split) & (saved.Representation == "ECFP") &
                           (saved.Threshold == threshold)].sort_values("TestOrder")
            assert np.array_equal(source.Index.to_numpy(int), te)
            original_model = SVC(C=1, kernel="precomputed", class_weight="balanced")
            original_model.fit(kernel[np.ix_(tr, tr)], y[tr])
            check_scores = original_model.decision_function(kernel[np.ix_(te, tr)])
            saved_scores = source.Score.to_numpy(float)
            error = float(np.max(np.abs(check_scores - saved_scores)))
            assert error < 1e-9
            checks.append({"Split": split, "Threshold": threshold, "MaxAbsoluteScoreDifference": error})
            reduced_model = SVC(C=1, kernel="precomputed", class_weight="balanced")
            reduced_model.fit(kernel[np.ix_(trk, trk)], y[trk])
            reduced_scores = reduced_model.decision_function(kernel[np.ix_(tek, trk)])
            for stage in STAGES:
                idx = te if stage == "original_all" else tek
                scores = saved_scores if stage == "original_all" else (saved_scores[mask] if stage == "frozen_reduced" else reduced_scores)
                weights = full_weights if stage == "original_all" else full_weights[:, mask]
                yy = y[idx]
                estimate = float(roc_auc_score(yy, scores))
                boot = tie_aware_weighted_auc(yy, scores, weights)
                estimates[stage, threshold], draws[stage, threshold] = estimate, boot
                lo, hi = percentile_ci(boot)
                auc_rows.append({"Split": split, "Stage": stage, "Threshold": threshold,
                                 "n_train": len(trk) if stage == "refit_reduced" else len(tr),
                                 "n_test": len(idx), "n_pos": int(yy.sum()), "n_neg": int(len(yy) - yy.sum()),
                                 "AUROC": estimate, "CI_low": lo, "CI_high": hi, "B_requested": B,
                                 "B_valid": int(np.isfinite(boot).sum())})
                for i, score in zip(idx, scores):
                    predictions.append({"Split": split, "Stage": stage, "Threshold": threshold,
                                        "Index": int(i), "CAS": study.df.iloc[i].CAS,
                                        "Y": int(y[i]), "Score": float(score)})
        gaps, gap_draws = {}, {}
        for stage in STAGES:
            observed = np.array([estimates[stage, severe] - estimates[stage, 100] for severe in (11, 1)])
            boot = np.column_stack([draws[stage, severe] - draws[stage, 100] for severe in (11, 1)])
            finite = np.isfinite(boot).all(axis=1)
            sim_lo, sim_hi = simultaneous_centered_ci(observed, boot[finite])
            for j, severe in enumerate((11, 1)):
                gaps[stage, severe] = observed[j]
                gap_draws[stage, severe] = boot[:, j]
                lo, hi = percentile_ci(boot[:, j])
                contrast_rows.append({"Split": split, "Stage": stage, "Severe": severe, "Gap": observed[j],
                                      "PercentileCI_low": lo, "PercentileCI_high": hi,
                                      "SimultaneousCI_low": sim_lo[j], "SimultaneousCI_high": sim_hi[j],
                                      "B_requested": B, "B_valid_joint_family": int(finite.sum())})
        for name, a, b in (("evaluation_only", "frozen_reduced", "original_all"),
                           ("refit_only", "refit_reduced", "frozen_reduced"),
                           ("total", "refit_reduced", "original_all")):
            observed = np.array([gaps[a, severe] - gaps[b, severe] for severe in (11, 1)])
            boot = np.column_stack([gap_draws[a, severe] - gap_draws[b, severe] for severe in (11, 1)])
            finite = np.isfinite(boot).all(axis=1)
            sim_lo, sim_hi = simultaneous_centered_ci(observed, boot[finite])
            for j, severe in enumerate((11, 1)):
                lo, hi = percentile_ci(boot[:, j])
                increments.append({"Split": split, "Component": name, "Severe": severe,
                                   "GapChange": observed[j], "PercentileCI_low": lo, "PercentileCI_high": hi,
                                   "SimultaneousCI_low": sim_lo[j], "SimultaneousCI_high": sim_hi[j],
                                   "B_requested": B, "B_valid_joint_family": int(finite.sum())})
        for severe in (11, 1):
            selected = {r["Component"]: r["GapChange"] for r in increments if r["Split"] == split and r["Severe"] == severe}
            assert np.isclose(selected["evaluation_only"] + selected["refit_only"], selected["total"], atol=1e-14)
    assert len(checks) == 9 and len(auc_rows) == 27 and len(contrast_rows) == len(increments) == 18
    prefix = "current_identifier_flagged_exclusion"
    for suffix, rows in (("baseline_checks", checks), ("cohort_counts", counts), ("predictions", predictions),
                         ("auroc", auc_rows), ("contrasts", contrast_rows), ("components", increments)):
        pd.DataFrame(rows).to_csv(OUT / f"{prefix}_{suffix}.csv", index=False)
    metadata = {"Status": "completed", "CohortName": prefix, "FlaggedN": 33,
                "FlagInterpretation": "current CID structure-connectivity flags, not 33 confirmed historical errors",
                "FlagSelectionRule": freeze["selection_rule"], "FlagFreezeUTC": freeze["frozen_utc"],
                "FlaggedCSV_SHA256": FLAG_HASH, "Model": "ECFP/Tanimoto SVC C=1, balanced class weights",
                "OriginalFitChecks": 9, "NewReducedFits": 9,
                "MaxBaselineScoreDifference": max(r["MaxAbsoluteScoreDifference"] for r in checks),
                "BootstrapDraws": B, "BootstrapSeedBase": SEED + 18000,
                "Bootstrap": "same scaffold draws on original test cohort for all endpoints and stages; reduced stages filter original draws",
                "Intervals": "percentile intervals and centered-bootstrap simultaneous intervals for the two severe-minus-broad contrasts within split/stage",
                "InferenceScope": "conditional and exploratory; deletion changes both membership and class weights; not causal error correction",
                "DataPolicy": "no replacement of SMILES, label, year or use flags; fixed official membership before explicit deletion",
                "Python": platform.python_version(), "Numpy": np.__version__, "Pandas": pd.__version__, "ScikitLearn": sklearn.__version__,
                "InputSHA256": {p: sha(ROOT / p) for p in ["data/raw/apistox.csv", "results/cache/kernel_ecfp.npy",
                     "results/cache/scaffolds.csv", "results/primary/02_test_predictions.csv"]}}
    (OUT / "CURRENT_IDENTIFIER_FLAGGED_EXCLUSION_METADATA.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(pd.DataFrame(counts).to_string(index=False))
    print(pd.DataFrame(auc_rows).pivot(index=["Split", "Stage"], columns="Threshold", values="AUROC").round(6).to_string())
    print(pd.DataFrame(contrast_rows).to_string(index=False))
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
