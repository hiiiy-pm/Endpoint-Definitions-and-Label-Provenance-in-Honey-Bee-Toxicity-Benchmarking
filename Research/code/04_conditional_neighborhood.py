from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.preprocessing import RobustScaler
from statsmodels.stats.multitest import multipletests

from _common import (
    RESULTS,
    SEED,
    build_or_load_cache,
    entropy_from_neighbor_labels,
    permute_within_groups,
)

OUT = RESULTS / "geometry"
OUT.mkdir(parents=True, exist_ok=True)

cache = build_or_load_cache(force=False)
study = cache["study"]
tier = study.tier.astype(np.int8, copy=False)
df = study.df
kernels = {k: np.asarray(v) for k, v in cache["kernels"].items()}

# Descriptor-space neighborhood uses the same robust scaling family as the primary
# descriptor model. Scaling is unsupervised and used only to define the frozen graph.
z = RobustScaler().fit_transform(np.asarray(cache["descriptors"]))
d2 = ((z[:, None, :] - z[None, :, :]) ** 2).sum(axis=2)
similarities = {**kernels, "Descriptors12": -d2}

source_insecticide = (
    df["source"].astype(str) + "|I=" + df["insecticide"].astype(str)
).to_numpy()
source_route_insecticide = (
    df["source"].astype(str)
    + "|"
    + df["toxicity_type"].astype(str)
    + "|I="
    + df["insecticide"].astype(str)
).to_numpy()

schemes = {
    "Global": None,
    "Source×Insecticide": source_insecticide,
    "Source×Route×Insecticide": source_route_insecticide,
}

P = 2000
BATCH = 200
K_VALUES = (5, 10, 20)
summary_rows: list[dict] = []
test_rows: list[dict] = []
null_payload: dict[str, np.ndarray] = {}


def permuted_batch(
    labels: np.ndarray,
    group_labels: np.ndarray | None,
    n: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Generate a compact batch of global or within-stratum permutations."""
    if group_labels is None:
        # argsort of random uniforms vectorizes independent row permutations and
        # avoids thousands of Python-level rng.permutation calls.
        order = np.argsort(rng.random((n, labels.size)), axis=1)
        return labels[order].astype(np.int8, copy=False)
    return np.stack(
        [permute_within_groups(labels, group_labels, rng) for _ in range(n)]
    ).astype(np.int8, copy=False)


def entropy_from_counts(counts: np.ndarray, k: int) -> np.ndarray:
    """Normalized four-class entropy from (..., 4) class counts."""
    p = counts.astype(np.float32, copy=False) / float(k)
    logp = np.zeros_like(p, dtype=np.float32)
    np.log(p, out=logp, where=p > 0)
    return -(p * logp).sum(axis=-1) / np.log(4.0)


rep_order = list(similarities)
scheme_order = list(schemes)

for rep_i, (rep, sim) in enumerate(similarities.items()):
    s = np.asarray(sim).copy()
    np.fill_diagonal(s, -np.inf)
    neighbors20 = np.argsort(-s, axis=1, kind="stable")[:, :20]

    observed_effects: dict[int, float] = {}
    for k in K_VALUES:
        neighbors = neighbors20[:, :k]
        h_obs = entropy_from_neighbor_labels(tier, neighbors)[0]
        same_obs = np.mean(tier[neighbors] == tier[:, None], axis=1)
        for t in range(4):
            mask = tier == t
            summary_rows.append(
                {
                    "Representation": rep,
                    "k": k,
                    "Tier": t,
                    "n": int(mask.sum()),
                    "MeanNormEntropy": float(h_obs[mask].mean()),
                    "MedianNormEntropy": float(np.median(h_obs[mask])),
                    "MeanSameTierNeighborFraction": float(same_obs[mask].mean()),
                }
            )
        observed_effects[k] = float(
            h_obs[tier == 2].mean() - h_obs[np.isin(tier, [1, 3])].mean()
        )

    # A single frozen permutation stream is shared across k within each
    # representation×conditioning scheme. This is statistically valid, improves
    # comparability across scales, and computes all three k values from one
    # 20-neighbor label tensor.
    for scheme_i, (scheme, group_labels) in enumerate(schemes.items()):
        rng = np.random.default_rng(SEED + rep_i * 100_000 + scheme_i * 1_000)
        null_parts: dict[int, list[np.ndarray]] = {k: [] for k in K_VALUES}

        for start in range(0, P, BATCH):
            n_batch = min(BATCH, P - start)
            perm = permuted_batch(tier, group_labels, n_batch, rng)
            neighbor_labels = perm[:, neighbors20]  # batch × molecule × 20

            # Cumulative four-class counts let k=5,10,20 reuse the same gather.
            cumulative_counts = np.stack(
                [np.cumsum(neighbor_labels == c, axis=2, dtype=np.int16) for c in range(4)],
                axis=3,
            )
            mask2 = perm == 2
            mask13 = (perm == 1) | (perm == 3)
            denom2 = mask2.sum(axis=1)
            denom13 = mask13.sum(axis=1)

            for k in K_VALUES:
                counts_k = cumulative_counts[:, :, k - 1, :]
                h_perm = entropy_from_counts(counts_k, k)
                effect = (
                    (h_perm * mask2).sum(axis=1) / denom2
                    - (h_perm * mask13).sum(axis=1) / denom13
                )
                null_parts[k].append(effect.astype(np.float32, copy=False))

        for k in K_VALUES:
            null = np.concatenate(null_parts[k]).astype(np.float64, copy=False)
            observed_effect = observed_effects[k]
            p_one_sided = float((1 + np.sum(null >= observed_effect)) / (len(null) + 1))
            test_rows.append(
                {
                    "Representation": rep,
                    "k": k,
                    "PermutationScheme": scheme,
                    "Comparison": "Tier2-vs-Tier1+Tier3",
                    "ObservedEffect": observed_effect,
                    "NullMean": float(null.mean()),
                    "NullSD": float(null.std(ddof=1)),
                    "NullQ025": float(np.quantile(null, 0.025)),
                    "NullQ975": float(np.quantile(null, 0.975)),
                    "P_one_sided": p_one_sided,
                    "Permutations": P,
                }
            )
            null_payload[f"{rep}_k{k}_{scheme}"] = null.astype(np.float32)

summary = pd.DataFrame(summary_rows)
summary["Representation"] = pd.Categorical(
    summary["Representation"], categories=rep_order, ordered=True
)
summary = summary.sort_values(["Representation", "k", "Tier"]).reset_index(drop=True)
summary["Representation"] = summary["Representation"].astype(str)
summary.to_csv(OUT / "01_neighborhood_entropy_summary.csv", index=False)

tests = pd.DataFrame(test_rows)
tests["Representation"] = pd.Categorical(
    tests["Representation"], categories=rep_order, ordered=True
)
tests["PermutationScheme"] = pd.Categorical(
    tests["PermutationScheme"], categories=scheme_order, ordered=True
)
tests = tests.sort_values(["Representation", "k", "PermutationScheme"]).reset_index(drop=True)
tests["BH_q_within_scheme"] = np.nan
for scheme, g in tests.groupby("PermutationScheme", observed=False):
    tests.loc[g.index, "BH_q_within_scheme"] = multipletests(
        g["P_one_sided"].to_numpy(), method="fdr_bh"
    )[1]
tests["Significant_q05"] = tests["BH_q_within_scheme"] < 0.05
tests["Representation"] = tests["Representation"].astype(str)
tests["PermutationScheme"] = tests["PermutationScheme"].astype(str)
tests.to_csv(OUT / "02_conditional_neighborhood_permutation.csv", index=False)
np.savez_compressed(OUT / "02_conditional_neighborhood_nulls.npz", **null_payload)

strata_rows = []
for scheme, group_labels in schemes.items():
    if group_labels is None:
        strata_rows.append({"PermutationScheme": scheme, "Stratum": "ALL", "n": len(tier)})
    else:
        for stratum, n in pd.Series(group_labels).value_counts().items():
            strata_rows.append(
                {"PermutationScheme": scheme, "Stratum": str(stratum), "n": int(n)}
            )
pd.DataFrame(strata_rows).to_csv(OUT / "03_permutation_stratum_sizes.csv", index=False)

print("Observed neighborhood summary at k=10:")
print(summary[summary["k"] == 10].round(4).to_string(index=False))
print("\nPermutation evidence by conditioning scheme:")
print(
    tests.groupby("PermutationScheme", observed=False)
    .agg(
        Significant=("Significant_q05", "sum"),
        Total=("Significant_q05", "size"),
        Max_q=("BH_q_within_scheme", "max"),
        Median_effect=("ObservedEffect", "median"),
    )
    .round(4)
    .to_string()
)
