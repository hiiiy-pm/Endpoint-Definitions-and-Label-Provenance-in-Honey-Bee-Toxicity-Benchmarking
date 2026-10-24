from __future__ import annotations

import math

import numpy as np
import pandas as pd
import patsy
from scipy.stats import norm
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from statsmodels.genmod.cov_struct import Exchangeable
from statsmodels.genmod.families import Binomial, Gaussian
from statsmodels.genmod.generalized_estimating_equations import GEE
from statsmodels.stats.multitest import multipletests

from _common import RESULTS, SEED, THRESHOLDS, build_or_load_cache

OUT = RESULTS / "applicability"
OUT.mkdir(parents=True, exist_ok=True)

cache = build_or_load_cache(force=False)
study = cache["study"]
kernel = np.asarray(cache["kernels"]["ECFP"])


def fit_svc_and_train_only_platt(
    tr: np.ndarray, te: np.ndarray, y: np.ndarray, seed: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ytr = y[tr]
    local = np.arange(len(tr))
    oof = np.empty(len(tr), dtype=float)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for fit_local, val_local in cv.split(local, ytr):
        fit_idx = tr[fit_local]
        val_idx = tr[val_local]
        model = SVC(C=1, kernel="precomputed", class_weight="balanced")
        model.fit(kernel[np.ix_(fit_idx, fit_idx)], y[fit_idx])
        oof[val_local] = model.decision_function(kernel[np.ix_(val_idx, fit_idx)])

    calibrator = LogisticRegression(C=1e6, solver="lbfgs", max_iter=5000)
    calibrator.fit(oof.reshape(-1, 1), ytr)

    model = SVC(C=1, kernel="precomputed", class_weight="balanced")
    model.fit(kernel[np.ix_(tr, tr)], ytr)
    score = model.decision_function(kernel[np.ix_(te, tr)])
    probability = calibrator.predict_proba(score.reshape(-1, 1))[:, 1]
    prediction = (score >= 0).astype(np.int8)
    return score, probability, prediction


stacked_parts = []
interaction_rows = []
slope_rows = []
joint_rows = []

for split_i, (split, (tr, te)) in enumerate(study.splits.items()):
    max_similarity = kernel[np.ix_(te, tr)].max(axis=1)
    novelty = 1.0 - max_similarity
    records = []
    for thr_i, thr in enumerate(THRESHOLDS):
        score, probability, prediction = fit_svc_and_train_only_platt(
            tr, te, study.y[thr], seed=SEED + split_i * 100 + thr_i
        )
        yy = study.y[thr][te]
        for j, idx in enumerate(te):
            records.append(
                {
                    "Split": split,
                    "Index": int(idx),
                    "Boundary": str(thr),
                    "MaxSimilarity": float(max_similarity[j]),
                    "Novelty": float(novelty[j]),
                    "Y": int(yy[j]),
                    "Score": float(score[j]),
                    "Probability": float(probability[j]),
                    "BrierLoss": float((yy[j] - probability[j]) ** 2),
                    "ClassificationError": int(prediction[j] != yy[j]),
                }
            )
    st = pd.DataFrame(records)
    st["NoveltyCentered"] = st["Novelty"] - st["Novelty"].mean()
    stacked_parts.append(st)

    for outcome, family in [("BrierLoss", Gaussian()), ("ClassificationError", Binomial())]:
        formula = (
            f'{outcome} ~ C(Boundary, Treatment(reference="100")) * NoveltyCentered'
        )
        ymat, xmat = patsy.dmatrices(formula, st, return_type="dataframe")
        fit = GEE(
            ymat,
            xmat,
            groups=st["Index"],
            family=family,
            cov_struct=Exchangeable(),
        ).fit()

        interaction_names = [
            name for name in fit.params.index if ":" in name and "NoveltyCentered" in name
        ]
        raw_p = [float(fit.pvalues[name]) for name in interaction_names]
        holm = multipletests(raw_p, method="holm")[1]
        for name, adjusted_p in zip(interaction_names, holm):
            boundary = "11" if "[T.11]" in name else "1"
            coef = float(fit.params[name])
            se = float(fit.bse[name])
            interaction_rows.append(
                {
                    "Split": split,
                    "Outcome": outcome,
                    "Contrast": f"{boundary} vs 100 novelty slope",
                    "SlopeDifference": coef,
                    "SE": se,
                    "CI_low": coef - 1.96 * se,
                    "CI_high": coef + 1.96 * se,
                    "RawP": float(fit.pvalues[name]),
                    "HolmP": float(adjusted_p),
                }
            )

        restriction = np.zeros((len(interaction_names), len(fit.params)))
        for row, name in enumerate(interaction_names):
            restriction[row, list(fit.params.index).index(name)] = 1.0
        wald = fit.wald_test(restriction, scalar=True)
        joint_rows.append(
            {
                "Split": split,
                "Outcome": outcome,
                "JointInteractionChi2": float(wald.statistic),
                "df": len(interaction_names),
                "P": float(wald.pvalue),
            }
        )

        if outcome == "BrierLoss":
            covariance = fit.cov_params().to_numpy()
            parameters = fit.params.to_numpy()
            names = list(fit.params.index)
            for boundary in ["100", "11", "1"]:
                contrast = np.zeros(len(parameters))
                contrast[names.index("NoveltyCentered")] = 1.0
                if boundary != "100":
                    interaction_name = [
                        n
                        for n in interaction_names
                        if f"[T.{boundary}]" in n
                    ][0]
                    contrast[names.index(interaction_name)] = 1.0
                slope = float(contrast @ parameters)
                variance = float(contrast @ covariance @ contrast)
                se = math.sqrt(max(variance, 0.0))
                slope_01 = 0.1 * slope
                se_01 = 0.1 * se
                z = slope_01 / se_01 if se_01 > 0 else np.nan
                slope_rows.append(
                    {
                        "Split": split,
                        "Boundary": int(boundary),
                        "Slope_per_0.1_novelty": slope_01,
                        "SE_per_0.1": se_01,
                        "CI_low": slope_01 - 1.96 * se_01,
                        "CI_high": slope_01 + 1.96 * se_01,
                        "P_slope": float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan,
                    }
                )

stacked = pd.concat(stacked_parts, ignore_index=True)
stacked.to_csv(OUT / "01_calibrated_predictions_and_novelty.csv", index=False)
interactions = pd.DataFrame(interaction_rows)
interactions.to_csv(OUT / "02_boundary_novelty_interactions.csv", index=False)
slopes = pd.DataFrame(slope_rows)
slopes.to_csv(OUT / "03_boundary_specific_brier_slopes.csv", index=False)
joints = pd.DataFrame(joint_rows)
joints.to_csv(OUT / "04_joint_interaction_tests.csv", index=False)

summary = (
    stacked.groupby(["Split", "Boundary"], as_index=False)
    .agg(
        n=("Index", "size"),
        MeanNovelty=("Novelty", "mean"),
        MeanBrier=("BrierLoss", "mean"),
        ErrorRate=("ClassificationError", "mean"),
    )
)
summary.to_csv(OUT / "05_applicability_summary.csv", index=False)

print("Boundary-specific Brier slopes per +0.1 novelty:")
print(slopes.round(5).to_string(index=False))
print("\nInteraction terms:")
print(interactions.round(5).to_string(index=False))
print("\nJoint interaction tests:")
print(joints.round(5).to_string(index=False))
