from __future__ import annotations

import numpy as np
import pandas as pd

from _common import RESULTS, SEED, build_or_load_cache, multinomial_weights, percentile_ci

OUT = RESULTS / "geometry"
OUT.mkdir(parents=True, exist_ok=True)

cache = build_or_load_cache(force=False)
study = cache["study"]
tier = study.tier
scaffolds = np.asarray(cache["scaffolds"], dtype=object)


def normalized_entropy(counts: np.ndarray) -> float:
    counts = np.asarray(counts, dtype=float)
    if counts.sum() <= 0:
        return np.nan
    p = counts[counts > 0] / counts.sum()
    return float(-(p * np.log(p)).sum() / np.log(4.0))


# ---------- full scaffold descriptors and leave-one-out molecular context ----------
scaffold_rows = []
molecule_rows = []
for scaffold in pd.unique(scaffolds):
    idx = np.flatnonzero(scaffolds == scaffold)
    counts = np.bincount(tier[idx], minlength=4)
    n = int(len(idx))
    full_entropy = normalized_entropy(counts)
    purity = float(counts.max() / n)
    n_tiers = int(np.sum(counts > 0))
    scaffold_rows.append(
        {
            "Scaffold": str(scaffold),
            "n": n,
            "T0": int(counts[0]),
            "T1": int(counts[1]),
            "T2": int(counts[2]),
            "T3": int(counts[3]),
            "FullPurity": purity,
            "FullNormEntropy": full_entropy,
            "n_tiers": n_tiers,
        }
    )
    for i in idx:
        loo_counts = counts.copy()
        loo_counts[tier[i]] -= 1
        molecule_rows.append(
            {
                "Index": int(i),
                "Tier": int(tier[i]),
                "Scaffold": str(scaffold),
                "Scaffold_n": n,
                "Scaffold_n_tiers": n_tiers,
                "LOO_Entropy": normalized_entropy(loo_counts) if n > 1 else np.nan,
                "LOO_Purity": float(loo_counts.max() / (n - 1)) if n > 1 else np.nan,
            }
        )

scaffold_df = pd.DataFrame(scaffold_rows).sort_values(["n", "Scaffold"], ascending=[False, True])
scaffold_df.to_csv(OUT / "04_scaffold_level_descriptors.csv", index=False)
molecule_df = pd.DataFrame(molecule_rows).sort_values("Index")
molecule_df.to_csv(OUT / "05_leave_one_out_scaffold_context.csv", index=False)

scaffold_tier = (
    molecule_df.dropna(subset=["LOO_Entropy"])
    .groupby(["Scaffold", "Scaffold_n", "Tier"], as_index=False)
    .agg(
        TierMolecule_n=("Index", "size"),
        MeanLOOEntropy=("LOO_Entropy", "mean"),
        MeanLOOPurity=("LOO_Purity", "mean"),
    )
)
scaffold_tier.to_csv(OUT / "06_scaffold_equal_weight_tier_context.csv", index=False)

# ---------- descriptive summaries ----------
summary_rows = []
for min_n in [3, 5]:
    eligible_mol = molecule_df[molecule_df["Scaffold_n"] >= min_n]
    eligible_st = scaffold_tier[scaffold_tier["Scaffold_n"] >= min_n]
    for t in range(4):
        gm = eligible_mol[eligible_mol["Tier"] == t]
        gs = eligible_st[eligible_st["Tier"] == t]
        summary_rows.append(
            {
                "MinScaffoldN": min_n,
                "Tier": t,
                "Molecule_n": int(len(gm)),
                "Scaffold_n": int(gs["Scaffold"].nunique()),
                "MoleculeWeightedMeanLOOEntropy": float(gm["LOO_Entropy"].mean()),
                "ScaffoldEqualMeanLOOEntropy": float(gs["MeanLOOEntropy"].mean()),
                "MoleculeWeightedMixedFraction": float(
                    (gm["Scaffold_n_tiers"] > 1).mean()
                ),
            }
        )
summary = pd.DataFrame(summary_rows)
summary.to_csv(OUT / "07_scaffold_context_summary.csv", index=False)

# ---------- scaffold-cluster bootstrap of leave-one-out context ----------
B = 5000
contrast_rows = []
npz_payload = {}

for min_i, min_n in enumerate([3, 5]):
    eligible_scaffolds = scaffold_df.loc[scaffold_df["n"] >= min_n, "Scaffold"].to_numpy()
    g = len(eligible_scaffolds)
    scaffold_to_row = {s: i for i, s in enumerate(eligible_scaffolds)}

    sum_entropy = np.zeros((g, 4), dtype=float)
    count_molecules = np.zeros((g, 4), dtype=float)
    mean_entropy = np.zeros((g, 4), dtype=float)
    present = np.zeros((g, 4), dtype=float)

    sub_mol = molecule_df[molecule_df["Scaffold"].isin(eligible_scaffolds)]
    for (scaffold, t), part in sub_mol.groupby(["Scaffold", "Tier"]):
        i = scaffold_to_row[scaffold]
        t = int(t)
        vals = part["LOO_Entropy"].dropna().to_numpy()
        if len(vals):
            sum_entropy[i, t] = vals.sum()
            count_molecules[i, t] = len(vals)
            mean_entropy[i, t] = vals.mean()
            present[i, t] = 1.0

    rng = np.random.default_rng(SEED + 5100 + min_i)
    weights = multinomial_weights(g, B, rng).astype(float)

    molecule_boot = np.full((B, 4), np.nan)
    scaffold_equal_boot = np.full((B, 4), np.nan)
    for t in range(4):
        numerator = weights @ sum_entropy[:, t]
        denominator = weights @ count_molecules[:, t]
        valid = denominator > 0
        molecule_boot[valid, t] = numerator[valid] / denominator[valid]

        numerator_eq = weights @ mean_entropy[:, t]
        denominator_eq = weights @ present[:, t]
        valid_eq = denominator_eq > 0
        scaffold_equal_boot[valid_eq, t] = numerator_eq[valid_eq] / denominator_eq[valid_eq]

    observed_molecule = np.array(
        [
            sub_mol.loc[sub_mol["Tier"] == t, "LOO_Entropy"].mean()
            for t in range(4)
        ]
    )
    eligible_st = scaffold_tier[scaffold_tier["Scaffold_n"] >= min_n]
    observed_equal = np.array(
        [
            eligible_st.loc[eligible_st["Tier"] == t, "MeanLOOEntropy"].mean()
            for t in range(4)
        ]
    )

    for weighting, boot_arr, observed_arr in [
        ("MoleculeWeighted", molecule_boot, observed_molecule),
        ("ScaffoldEqualWeighted", scaffold_equal_boot, observed_equal),
    ]:
        for other in [0, 1, 3]:
            d = boot_arr[:, 2] - boot_arr[:, other]
            d = d[np.isfinite(d)]
            low, high = percentile_ci(d)
            contrast_rows.append(
                {
                    "MinScaffoldN": min_n,
                    "Analysis": weighting,
                    "Comparison": f"T2-T{other}",
                    "ObservedDifference": float(observed_arr[2] - observed_arr[other]),
                    "BootstrapMean": float(d.mean()),
                    "CI_low": low,
                    "CI_high": high,
                    "ProbDifference_gt0": float(np.mean(d > 0)),
                    "B_valid": int(len(d)),
                    "ScaffoldsResampled": g,
                }
            )
            npz_payload[f"min{min_n}_{weighting}_T2_T{other}"] = d.astype(np.float32)

    # Paired-scaffold sensitivity: only scaffolds that contain both Tier2 and comparator.
    wide = eligible_st.pivot(index="Scaffold", columns="Tier", values="MeanLOOEntropy")
    for other in [0, 1, 3]:
        paired = wide[[2, other]].dropna()
        diff = (paired[2] - paired[other]).to_numpy()
        if not len(diff):
            continue
        rng_pair = np.random.default_rng(SEED + 6100 + min_i * 10 + other)
        w_pair = multinomial_weights(len(diff), B, rng_pair).astype(float)
        boot_diff = (w_pair @ diff) / w_pair.sum(axis=1)
        low, high = percentile_ci(boot_diff)
        contrast_rows.append(
            {
                "MinScaffoldN": min_n,
                "Analysis": "PairedScaffoldsContainingBothTiers",
                "Comparison": f"T2-T{other}",
                "ObservedDifference": float(diff.mean()),
                "BootstrapMean": float(boot_diff.mean()),
                "CI_low": low,
                "CI_high": high,
                "ProbDifference_gt0": float(np.mean(boot_diff > 0)),
                "B_valid": int(len(boot_diff)),
                "ScaffoldsResampled": int(len(diff)),
            }
        )
        npz_payload[f"min{min_n}_Paired_T2_T{other}"] = boot_diff.astype(np.float32)

contrasts = pd.DataFrame(contrast_rows)
contrasts.to_csv(OUT / "08_scaffold_sensitivity_bootstrap.csv", index=False)
np.savez_compressed(OUT / "08_scaffold_sensitivity_bootstrap_distributions.npz", **npz_payload)

print("Leave-one-out scaffold context summary:")
print(summary.round(4).to_string(index=False))
print("\nScaffold-aware contrasts:")
print(contrasts.round(4).to_string(index=False))
