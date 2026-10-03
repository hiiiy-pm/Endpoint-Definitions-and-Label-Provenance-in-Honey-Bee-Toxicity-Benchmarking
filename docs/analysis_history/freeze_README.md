# ApisTox Final Paper Package v1

## Purpose

This package freezes the final, paper-ready analysis after the exploratory, mechanism-validation, and final statistical-seal rounds. It is designed so that the manuscript can be written without mixing superseded claims from earlier rounds.

## Final research positioning

**Recommended working title**

> Regulatory-Boundary Heterogeneity in Machine-Learning Prediction of Honey Bee Acute Toxicity: Chemical-Space Ambiguity and Out-of-Distribution Generalization

**Chinese working title**

> 蜜蜂急性毒性机器学习预测中的监管边界异质性：化学空间模糊性与分布外泛化

The paper is **not** positioned as a new-model/SOTA paper. Its main contribution is a scientific finding: molecular learnability differs across regulatory toxicity boundaries, with the 100 μg/bee boundary consistently weaker than the 11 and 1 μg/bee severe-toxicity boundaries, accompanied by local/scaffold chemical-space ambiguity.

## Frozen primary findings

1. **Regulatory-boundary heterogeneity — strongly supported.** The 100 μg/bee endpoint is consistently less learnable than the 11 and 1 μg/bee endpoints across 5 representations × 3 official splits (15/15 settings).
2. **Chemical-space ambiguity — supported mechanism.** Tier2 shows persistent excess local neighborhood mixing across representations and neighborhood sizes; scaffold-level analyses support higher Tier2 ambiguity relative to Tier1, with sensitivity limits versus Tier3 at stricter scaffold-size filtering.
3. **Chemical novelty — secondary applicability-domain result.** Novelty increases prediction error under OOD splits, especially for the 100 μg endpoint; however, formal Boundary × Novelty interaction tests are not significant and must not be presented as an established mechanism.

## Claims that are NOT supported

- “1 μg always outperforms 11 μg.”
- “Tier2 is always the single most ambiguous tier.”
- “Boundary × chemical-novelty interaction is statistically significant.”
- “Time split conclusively confirms all boundary differences.”
- “A compact physicochemical signature transfers robustly through time.”

See `docs/FINAL_CLAIMS.md` for exact wording.

## Package structure

- `baseline_round2/` — frozen second-round baseline, raw ApisTox CSV, reconstructed official splits, historical comparison outputs.
- `code/` — portable final-analysis scripts using relative paths.
- `results/stage2/` — representation robustness, novelty strata, scaffold overlap, negative control, descriptor effects.
- `results/mechanism/` — 5000 bootstraps, neighborhood entropy, scaffold analyses, sparse descriptor experiments.
- `results/final_seal/` — GEE interaction tests, multiplicity correction, robustness matrix, final evidence matrix.
- `docs/` — paper-writing master notes, final claim boundaries, result-file map, methods checklist.

## Reproduction

### Conda

```bash
conda env create -f environment.yml
conda activate apistox-final-paper
python code/run_all.py
```

### pip

```bash
python -m venv .venv
# Windows: .venv\\Scripts\\activate
# Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
python code/run_all.py
```

The scripts use fixed seed `20260829` where stochastic procedures are involved. The historical round-2 baseline is frozen and is not rerun by `run_all.py`.

## Main data definition

The source file contains 1035 molecules. The four derived tiers are deterministically reconstructed from `label` and `ppdb_level`:

- Tier0: >100 μg/bee, n=177
- Tier1: 11–100 μg/bee, n=562
- Tier2: 1–11 μg/bee, n=125
- Tier3: ≤1 μg/bee, n=171

Three cumulative endpoints are then defined:

- `Y100 = 1(Tier >= 1)` → LD50 ≤100 μg/bee
- `Y11  = 1(Tier >= 2)` → LD50 ≤11 μg/bee
- `Y1   = 1(Tier >= 3)` → LD50 ≤1 μg/bee

Official reconstructed splits are 828 train / 207 test for Random, MaxMin, and Time.

## Integrity note

No continuous LD50 values are contained in the frozen source package, so threshold-proximity analysis was deliberately not performed. Continuous LD50 values must not be reverse-engineered from the categorical labels.
