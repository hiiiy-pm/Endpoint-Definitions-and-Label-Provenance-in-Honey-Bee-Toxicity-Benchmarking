from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests

from _common import (
    RAW_CSV,
    RESULTS,
    SEED,
    THRESHOLDS,
    build_or_load_cache,
    cramer_v,
    sha256_file,
)

OUT = RESULTS / "decomposition"
OUT.mkdir(parents=True, exist_ok=True)

cache = build_or_load_cache(force=False)
study = cache["study"]
df = study.df.copy()
df["Tier"] = study.tier

# ---------- hard validation ----------
validation = {
    "seed": SEED,
    "n_rows": int(len(df)),
    "n_unique_smiles": int(df["SMILES"].nunique()),
    "tier_counts": {str(k): int(v) for k, v in df["Tier"].value_counts().sort_index().items()},
    "source_counts": {str(k): int(v) for k, v in df["source"].value_counts().items()},
    "nested_endpoint_invariant": bool(
        np.all(study.y[1] <= study.y[11]) and np.all(study.y[11] <= study.y[100])
    ),
    "source_csv_sha256": sha256_file(RAW_CSV),
    "cache_validation": cache["validation"],
}
for split, (tr, te) in study.splits.items():
    validation[f"split_{split}"] = {
        "train_n": int(len(tr)),
        "test_n": int(len(te)),
        "overlap_n": int(len(set(tr).intersection(te))),
        "train_unique_n": int(len(np.unique(tr))),
        "test_unique_n": int(len(np.unique(te))),
        "union_n": int(len(set(tr).union(te))),
    }
(RESULTS / "final").mkdir(parents=True, exist_ok=True)
(RESULTS / "final" / "00_data_validation.json").write_text(
    json.dumps(validation, indent=2, ensure_ascii=False), encoding="utf-8"
)

# ---------- endpoint and split class counts ----------
rows = []
for split, (tr, te) in study.splits.items():
    for subset, idx in [("Train", tr), ("Test", te)]:
        for thr in THRESHOLDS:
            yy = study.y[thr][idx]
            rows.append(
                {
                    "Split": split,
                    "Subset": subset,
                    "Threshold": thr,
                    "n": len(idx),
                    "Positive_n": int(yy.sum()),
                    "Negative_n": int(len(yy) - yy.sum()),
                    "Positive_fraction": float(yy.mean()),
                }
            )
pd.DataFrame(rows).to_csv(OUT / "00_split_endpoint_counts.csv", index=False)

# ---------- tier composition ----------
# Categorical variables are represented by mutually exclusive levels; agrochemical flags
# are separate, non-mutually-exclusive binary indicators.
long_rows = []
for tier, g in df.groupby("Tier", sort=True):
    for col in ["source", "toxicity_type"]:
        counts = g[col].value_counts(dropna=False)
        for level, n in counts.items():
            long_rows.append(
                {
                    "Tier": int(tier),
                    "Variable": col,
                    "Level": str(level),
                    "n": int(n),
                    "Tier_n": int(len(g)),
                    "Proportion": float(n / len(g)),
                }
            )
    for col in ["herbicide", "fungicide", "insecticide", "other_agrochemical"]:
        n = int(g[col].sum())
        long_rows.append(
            {
                "Tier": int(tier),
                "Variable": "agrochemical_flag",
                "Level": col,
                "n": n,
                "Tier_n": int(len(g)),
                "Proportion": float(n / len(g)),
            }
        )
comp_long = pd.DataFrame(long_rows)
comp_long.to_csv(OUT / "01_tier_composition_long.csv", index=False)
comp_wide = comp_long.pivot_table(
    index="Tier", columns=["Variable", "Level"], values="Proportion", aggfunc="first"
).sort_index()
comp_wide.columns = [f"{a}:{b}" for a, b in comp_wide.columns]
comp_wide.reset_index().to_csv(OUT / "02_tier_composition_wide.csv", index=False)

# ---------- association tests ----------
assoc_rows = []
for col in [
    "source",
    "toxicity_type",
    "herbicide",
    "fungicide",
    "insecticide",
    "other_agrochemical",
]:
    tab = pd.crosstab(df["Tier"], df[col])
    chi2, p, dof, _ = chi2_contingency(tab)
    assoc_rows.append(
        {
            "Variable": col,
            "Chi2": float(chi2),
            "df": int(dof),
            "P": float(p),
            "CramerV": cramer_v(tab),
        }
    )
assoc = pd.DataFrame(assoc_rows)
assoc["BH_q"] = multipletests(assoc["P"], method="fdr_bh")[1]
assoc.to_csv(OUT / "03_tier_composition_association_tests.csv", index=False)

# ---------- human-readable audit summary ----------
lines = [
    "# Revised data and composition audit",
    "",
    f"- Frozen molecules: {len(df):,}",
    f"- Tier counts: {validation['tier_counts']}",
    f"- Source counts: {validation['source_counts']}",
    f"- Nested endpoint invariant: {validation['nested_endpoint_invariant']}",
    "",
    "## Composition proportions by tier",
    "",
    comp_wide.reset_index().to_markdown(index=False, floatfmt=".3f"),
    "",
    "## Tier-composition association tests",
    "",
    assoc.to_markdown(index=False, floatfmt=".6g"),
]
(OUT / "DATA_COMPOSITION_SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")

print(json.dumps(validation, indent=2, ensure_ascii=False))
print("\nTier composition (proportions):")
print(comp_wide.round(3).to_string())
print("\nAssociation tests:")
print(assoc.round(6).to_string(index=False))
