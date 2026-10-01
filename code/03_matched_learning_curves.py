from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.svm import SVC

from _common import RESULTS, SEED, THRESHOLDS, build_or_load_cache, percentile_ci

OUT = RESULTS / "learning_curves"
OUT.mkdir(parents=True, exist_ok=True)

cache = build_or_load_cache(force=False)
study = cache["study"]
kernel = np.asarray(cache["kernels"]["ECFP"])

PER_CLASS_SIZES = [25, 50, 75, 100, 125]
REPEATS = 200
rows = []

for split_i, (split, (tr, te)) in enumerate(study.splits.items()):
    for threshold_i, thr in enumerate(THRESHOLDS):
        y = study.y[thr]
        pos = tr[y[tr] == 1]
        neg = tr[y[tr] == 0]
        for size_i, n_per_class in enumerate(PER_CLASS_SIZES):
            if n_per_class > min(len(pos), len(neg)):
                continue
            rng = np.random.default_rng(
                SEED + split_i * 100_000 + threshold_i * 10_000 + size_i * 1_000
            )
            for repeat in range(REPEATS):
                selected = np.concatenate(
                    [
                        rng.choice(pos, n_per_class, replace=False),
                        rng.choice(neg, n_per_class, replace=False),
                    ]
                )
                rng.shuffle(selected)
                model = SVC(
                    C=1,
                    kernel="precomputed",
                    class_weight="balanced",
                    cache_size=1024,
                )
                model.fit(kernel[np.ix_(selected, selected)], y[selected])
                score = model.decision_function(kernel[np.ix_(te, selected)])
                rows.append(
                    {
                        "Split": split,
                        "Threshold": thr,
                        "n_per_class": n_per_class,
                        "Train_n": int(2 * n_per_class),
                        "Repeat": repeat,
                        "AUROC": float(roc_auc_score(y[te], score)),
                        "AUPRC": float(average_precision_score(y[te], score)),
                    }
                )

long = pd.DataFrame(rows)
long.to_csv(OUT / "01_matched_learning_curves_long.csv", index=False)

summary_rows = []
for keys, g in long.groupby(["Split", "Threshold", "n_per_class"], sort=False):
    split, thr, n_per_class = keys
    auc_low, auc_high = percentile_ci(g["AUROC"].to_numpy())
    ap_low, ap_high = percentile_ci(g["AUPRC"].to_numpy())
    summary_rows.append(
        {
            "Split": split,
            "Threshold": int(thr),
            "n_per_class": int(n_per_class),
            "Train_n": int(2 * n_per_class),
            "Repeats": int(len(g)),
            "MeanAUROC": float(g["AUROC"].mean()),
            "SDAUROC": float(g["AUROC"].std(ddof=1)),
            "AUROC_Q025": auc_low,
            "AUROC_Q975": auc_high,
            "MeanAUPRC": float(g["AUPRC"].mean()),
            "AUPRC_Q025": ap_low,
            "AUPRC_Q975": ap_high,
        }
    )
summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT / "02_matched_learning_curves_summary.csv", index=False)

# Descriptive endpoint contrasts over the repeated matched-sample fits.
contrast_rows = []
for (split, n_per_class), g in long.groupby(["Split", "n_per_class"]):
    pivot = g.pivot(index="Repeat", columns="Threshold", values="AUROC")
    for a, b in [(11, 100), (1, 100), (1, 11)]:
        d = (pivot[a] - pivot[b]).dropna().to_numpy()
        low, high = percentile_ci(d)
        contrast_rows.append(
            {
                "Split": split,
                "n_per_class": int(n_per_class),
                "Comparison": f"{a}-{b}",
                "MeanDeltaAUROC": float(d.mean()),
                "Q025": low,
                "Q975": high,
                "ProbDelta_gt0": float(np.mean(d > 0)),
                "Repeats": int(len(d)),
            }
        )
contrasts = pd.DataFrame(contrast_rows)
contrasts.to_csv(OUT / "03_matched_learning_curve_contrasts.csv", index=False)

print("Matched learning-curve mean AUROC:")
print(
    summary.pivot_table(
        index=["Split", "n_per_class"], columns="Threshold", values="MeanAUROC"
    ).round(3).to_string()
)
print("\nAt n=125 per class:")
print(contrasts[contrasts["n_per_class"] == 125].round(4).to_string(index=False))
