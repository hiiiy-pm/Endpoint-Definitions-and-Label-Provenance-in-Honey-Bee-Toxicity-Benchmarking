from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import RobustScaler
from sklearn.svm import SVC

from _common import (
    REPRESENTATIONS,
    RESULTS,
    SEED,
    THRESHOLDS,
    build_or_load_cache,
    multinomial_weights,
    percentile_ci,
    simultaneous_centered_ci,
    tie_aware_weighted_auc,
)

OUT = RESULTS / "primary"
OUT.mkdir(parents=True, exist_ok=True)

cache = build_or_load_cache(force=False)
study = cache["study"]
kernels = cache["kernels"]
descriptors = np.asarray(cache["descriptors"])
scaffolds = np.asarray(cache["scaffolds"], dtype=object)

metric_rows = []
pred_rows = []
pred_store: dict[tuple[str, str, int], np.ndarray] = {}

for split, (tr, te) in study.splits.items():
    for rep, kernel in kernels.items():
        ktr = np.asarray(kernel[np.ix_(tr, tr)])
        kte = np.asarray(kernel[np.ix_(te, tr)])
        for thr in THRESHOLDS:
            y = study.y[thr]
            model = SVC(C=1, kernel="precomputed", class_weight="balanced")
            model.fit(ktr, y[tr])
            score = model.decision_function(kte)
            pred_store[(split, rep, thr)] = score
            metric_rows.append(
                {
                    "Split": split,
                    "Representation": rep,
                    "Threshold": thr,
                    "AUROC": float(roc_auc_score(y[te], score)),
                    "AUPRC": float(average_precision_score(y[te], score)),
                    "n_test": int(len(te)),
                    "n_pos": int(y[te].sum()),
                }
            )
            for local_i, global_i in enumerate(te):
                pred_rows.append(
                    {
                        "Split": split,
                        "Representation": rep,
                        "Threshold": thr,
                        "TestOrder": local_i,
                        "Index": int(global_i),
                        "CID": int(study.df.iloc[global_i]["CID"]),
                        "Scaffold": str(scaffolds[global_i]),
                        "Y": int(y[global_i]),
                        "Score": float(score[local_i]),
                    }
                )

    scaler = RobustScaler().fit(descriptors[tr])
    xtr = scaler.transform(descriptors[tr])
    xte = scaler.transform(descriptors[te])
    for thr in THRESHOLDS:
        y = study.y[thr]
        model = SVC(C=1, kernel="rbf", gamma="scale", class_weight="balanced")
        model.fit(xtr, y[tr])
        score = model.decision_function(xte)
        rep = "Descriptors12"
        pred_store[(split, rep, thr)] = score
        metric_rows.append(
            {
                "Split": split,
                "Representation": rep,
                "Threshold": thr,
                "AUROC": float(roc_auc_score(y[te], score)),
                "AUPRC": float(average_precision_score(y[te], score)),
                "n_test": int(len(te)),
                "n_pos": int(y[te].sum()),
            }
        )
        for local_i, global_i in enumerate(te):
            pred_rows.append(
                {
                    "Split": split,
                    "Representation": rep,
                    "Threshold": thr,
                    "TestOrder": local_i,
                    "Index": int(global_i),
                    "CID": int(study.df.iloc[global_i]["CID"]),
                    "Scaffold": str(scaffolds[global_i]),
                    "Y": int(y[global_i]),
                    "Score": float(score[local_i]),
                }
            )

metrics = pd.DataFrame(metric_rows)
metrics.to_csv(OUT / "01_representation_performance.csv", index=False)
predictions = pd.DataFrame(pred_rows)
predictions.to_csv(OUT / "02_test_predictions.csv", index=False)

# ---------- Paired molecule/scaffold-cluster bootstrap for every representation ----------
# ECFP is the prespecified primary structural probe. The identical resampling weights are
# reused across representations so that the confirmatory analyses are paired with the
# primary one and ECFP estimates remain numerically unchanged.
B = 5000
contrasts = [(11, 100), (1, 100), (1, 11)]
primary_contrasts = [(11, 100), (1, 100)]
PRIMARY_REPRESENTATION = "ECFP"
boot_rows = []
npz_payload = {}
npz_payload_all = {}

for split_i, (split, (tr, te)) in enumerate(study.splits.items()):
    n = len(te)
    rng = np.random.default_rng(SEED + 101 * (split_i + 1))
    molecule_weights = multinomial_weights(n, B, rng)

    test_scaffolds = scaffolds[te]
    unique_scaffolds, inverse = np.unique(test_scaffolds, return_inverse=True)
    scaffold_counts = multinomial_weights(len(unique_scaffolds), B, rng)
    scaffold_weights = scaffold_counts[:, inverse]

    for mode, weights in [
        ("Molecule", molecule_weights),
        ("ScaffoldCluster", scaffold_weights),
    ]:
        for rep in REPRESENTATIONS:
            auc_boot = {}
            auc_obs = {}
            for thr in THRESHOLDS:
                yy = study.y[thr][te]
                ss = pred_store[(split, rep, thr)]
                auc_obs[thr] = float(roc_auc_score(yy, ss))
                auc_boot[thr] = tie_aware_weighted_auc(yy, ss, weights)

            primary_obs = np.array([auc_obs[a] - auc_obs[b] for a, b in primary_contrasts])
            primary_boot = np.column_stack(
                [auc_boot[a] - auc_boot[b] for a, b in primary_contrasts]
            )
            finite_primary = np.all(np.isfinite(primary_boot), axis=1)
            sim_low, sim_high = simultaneous_centered_ci(
                primary_obs, primary_boot[finite_primary], alpha=0.05
            )

            for a, b in contrasts:
                values = auc_boot[a] - auc_boot[b]
                values = values[np.isfinite(values)]
                obs = auc_obs[a] - auc_obs[b]
                low, high = percentile_ci(values)
                is_primary = (a, b) in primary_contrasts
                if is_primary:
                    j = primary_contrasts.index((a, b))
                    family_low, family_high = float(sim_low[j]), float(sim_high[j])
                else:
                    family_low, family_high = np.nan, np.nan
                boot_rows.append(
                    {
                        "Representation": rep,
                        "Split": split,
                        "Bootstrap": mode,
                        "Comparison": f"{a}-{b}",
                        "PrimaryComparison": is_primary,
                        "B_valid": int(len(values)),
                        "ObservedDeltaAUROC": float(obs),
                        "BootstrapMean": float(values.mean()),
                        "PercentileCI_low": low,
                        "PercentileCI_high": high,
                        "SimultaneousCI_low": family_low,
                        "SimultaneousCI_high": family_high,
                        "ProbDelta_gt0": float(np.mean(values > 0)),
                    }
                )
                npz_payload_all[f"{rep}_{split}_{mode}_{a}_{b}"] = values.astype(np.float32)
                if rep == PRIMARY_REPRESENTATION:
                    npz_payload[f"{split}_{mode}_{a}_{b}"] = values.astype(np.float32)

bootstrap_all = pd.DataFrame(boot_rows)
bootstrap_all.to_csv(OUT / "05_representation_paired_bootstrap.csv", index=False)
np.savez_compressed(
    OUT / "05_representation_paired_bootstrap_distributions.npz", **npz_payload_all
)

bootstrap = (
    bootstrap_all[bootstrap_all["Representation"] == PRIMARY_REPRESENTATION]
    .drop(columns=["Representation"])
    .reset_index(drop=True)
)
bootstrap.to_csv(OUT / "03_ecfp_paired_bootstrap.csv", index=False)
np.savez_compressed(OUT / "03_ecfp_paired_bootstrap_distributions.npz", **npz_payload)

# Descriptive cross-representation endpoint ordering count.
ordering = []
for (split, rep), g in metrics.groupby(["Split", "Representation"], sort=False):
    d = g.set_index("Threshold")["AUROC"]
    ordering.append(
        {
            "Split": split,
            "Representation": rep,
            "AUC100_lt_AUC11": bool(d.loc[100] < d.loc[11]),
            "AUC100_lt_AUC1": bool(d.loc[100] < d.loc[1]),
            "Both": bool(d.loc[100] < d.loc[11] and d.loc[100] < d.loc[1]),
        }
    )
ordering_df = pd.DataFrame(ordering)
ordering_df.to_csv(OUT / "04_endpoint_ordering_robustness.csv", index=False)

# Confirmatory cross-representation support for the two primary contrasts.
primary_rows = bootstrap_all[bootstrap_all["PrimaryComparison"]].copy()
primary_rows["CI_excludes_zero"] = primary_rows["SimultaneousCI_low"] > 0
support = (
    primary_rows.groupby(["Bootstrap", "Representation", "Split"], sort=False)["CI_excludes_zero"]
    .sum()
    .reset_index(name="n_primary_contrasts_supported")
)
support.to_csv(OUT / "06_cross_representation_primary_support.csv", index=False)

scaffold_support = primary_rows[primary_rows["Bootstrap"] == "ScaffoldCluster"]
summary = {
    "settings_total": int(len(ordering_df)),
    "settings_with_auc100_below_both_severe_endpoints": int(ordering_df["Both"].sum()),
    "bootstrap_repetitions": B,
    "primary_inference": "effect sizes with percentile and simultaneous centered-bootstrap confidence intervals",
    "primary_representation": PRIMARY_REPRESENTATION,
    "confirmatory_representations": [r for r in REPRESENTATIONS if r != PRIMARY_REPRESENTATION],
    "scaffold_cluster_primary_contrasts_total": int(len(scaffold_support)),
    "scaffold_cluster_primary_contrasts_ci_above_zero": int(
        (scaffold_support["SimultaneousCI_low"] > 0).sum()
    ),
}
(OUT / "PRIMARY_ANALYSIS_METADATA.json").write_text(
    json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
)

print(metrics.pivot_table(index=["Split", "Representation"], columns="Threshold", values="AUROC").round(3))
print("\nPrimary ECFP bootstrap:")
print(
    bootstrap[bootstrap["PrimaryComparison"]][
        [
            "Split",
            "Bootstrap",
            "Comparison",
            "ObservedDeltaAUROC",
            "PercentileCI_low",
            "PercentileCI_high",
            "SimultaneousCI_low",
            "SimultaneousCI_high",
            "ProbDelta_gt0",
        ]
    ].round(4).to_string(index=False)
)

print("\nScaffold-cluster primary contrasts, all representations:")
print(
    scaffold_support[
        [
            "Representation",
            "Split",
            "Comparison",
            "ObservedDeltaAUROC",
            "SimultaneousCI_low",
            "SimultaneousCI_high",
        ]
    ].round(4).to_string(index=False)
)
