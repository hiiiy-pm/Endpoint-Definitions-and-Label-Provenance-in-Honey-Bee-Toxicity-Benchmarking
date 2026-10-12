"""Figure: all four label-crossover cells and their paired full-revision effect.

The primary keep ECFP scenario is selected before inspecting effect sizes.
Every split and both severe endpoints are retained. The left panel gives point
estimates; the right panel supplies the paired uncertainty for the question.
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
QA = ROOT / "supplementary_audits/statistics/figure_qa"
QA.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(QA))
from audit_panel_alignment import require_matplotlib_panel_alignment
sys.path.insert(0, str(ROOT / "code"))
from _common import multinomial_weights, percentile_ci, tie_aware_weighted_auc

INPUT = ROOT / "results/label_audit/07_crossover_gap.csv"
PRED = ROOT / "results/label_audit/09_crossover_predictions.csv.gz"
SPLITS = ("Random", "MaxMin", "Time")
SEVERE = (11, 1)
CELLS = (("original", "original"), ("original", "corrected"),
         ("corrected", "original"), ("corrected", "corrected"))


def main():
    gap = pd.read_csv(INPUT)
    gap = gap[(gap.Scenario == "primary") & (gap.Unresolved == "keep") & (gap.Representation == "ECFP")]
    selected = []
    for split in SPLITS:
        for severe in SEVERE:
            rows = gap[(gap.Split == split) & (gap.Severe == severe)]
            for a, b in CELLS:
                row = rows[(rows.TrainLabels == a) & (rows.EvalLabels == b)]
                assert len(row) == 1
                selected.append({**row.iloc[0].to_dict(), "Panel": "a", "Quantity": "cell_gap"})
            change = rows[(rows.TrainLabels == "change:corrected-vs-original") &
                          (rows.EvalLabels == "change:corrected-vs-original")]
            assert len(change) == 1
            selected.append({**change.iloc[0].to_dict(), "Panel": "b", "Quantity": "paired_full_change"})
    source = pd.DataFrame(selected)
    assert len(source) == 30 and (source.n_train == 828).all() and (source.n_test == 207).all()
    assert (source.B_requested == 2000).all() and (source.B_valid == 2000).all()
    for _, rows in source.groupby("Split"):
        assert rows.train_index_sha256.nunique() == rows.test_index_sha256.nunique() == 1

    # Recompute from compound scores, including the exact published change draws.
    predictions = pd.read_csv(PRED)
    predictions = predictions[(predictions.Scenario == "primary") & (predictions.Unresolved == "keep") &
                              (predictions.Representation == "ECFP")]
    scaffolds = pd.read_csv(ROOT / "results/cache/scaffolds.csv").set_index("Index").Scaffold
    check_rows = []
    for split_index, split in enumerate(SPLITS):
        est, draws = {}, {}
        reference = predictions[(predictions.Split == split) & (predictions.Threshold == 100)]
        indices = reference.Index.to_numpy(int)
        uniq, inverse = np.unique(scaffolds.loc[indices].to_numpy(str), return_inverse=True)
        rng = np.random.default_rng(20260829 + 13000 + 97 * split_index)
        weights = multinomial_weights(len(uniq), 2000, rng)[:, inverse]
        for threshold in (100, 11, 1):
            rows = predictions[(predictions.Split == split) & (predictions.Threshold == threshold)]
            assert np.array_equal(rows.Index.to_numpy(int), indices)
            for train, evaluate in CELLS:
                scores = rows.OriginalTrainScore.to_numpy(float) if train == "original" else rows.CorrectedTrainScore.to_numpy(float)
                labels = rows.OriginalLabel.to_numpy(int) if evaluate == "original" else rows.CorrectedLabel.to_numpy(int)
                est[threshold, train, evaluate] = roc_auc_score(labels, scores)
                draws[threshold, train, evaluate] = tie_aware_weighted_auc(labels, scores, weights)
        for severe in SEVERE:
            for train, evaluate in CELLS:
                observed = est[severe, train, evaluate] - est[100, train, evaluate]
                row = source[(source.Split == split) & (source.Severe == severe) &
                             (source.TrainLabels == train) & (source.EvalLabels == evaluate)].iloc[0]
                difference = abs(observed - row.Gap)
                assert difference < 1e-12
                check_rows.append({"Split": split, "Severe": severe, "Quantity": f"{train}/{evaluate}",
                                   "MaxAbsoluteDifference": difference})
            diff_draws = ((draws[severe, "corrected", "corrected"] - draws[100, "corrected", "corrected"])
                          - (draws[severe, "original", "original"] - draws[100, "original", "original"]))
            lo, hi = percentile_ci(diff_draws)
            observed = ((est[severe, "corrected", "corrected"] - est[100, "corrected", "corrected"])
                        - (est[severe, "original", "original"] - est[100, "original", "original"]))
            row = source[(source.Split == split) & (source.Severe == severe) & (source.Quantity == "paired_full_change")].iloc[0]
            difference = max(abs(observed - row.Gap), abs(lo - row.CI_low), abs(hi - row.CI_high))
            assert difference < 1e-12
            check_rows.append({"Split": split, "Severe": severe, "Quantity": "paired_full_change_and_CI",
                               "MaxAbsoluteDifference": difference})
    source.to_csv(ROOT / "source_data/figure_label_crossover.csv", index=False)
    pd.DataFrame(check_rows).to_csv(QA / "numeric_verification.csv", index=False)

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "axes.unicode_minus": False,
                         "axes.linewidth": .65})
    fig, (left, right) = plt.subplots(1, 2, figsize=(7.16, 3.8), gridspec_kw={"width_ratios": [1, 1.1]})
    fig.subplots_adjust(left=.21, right=.98, bottom=.23, top=.80, wspace=.20)
    row_keys = [(split, severe) for split in SPLITS for severe in SEVERE]
    matrix = np.array([[source[(source.Split == split) & (source.Severe == severe) &
                                      (source.TrainLabels == a) & (source.EvalLabels == b)].iloc[0].Gap
                        for a, b in CELLS] for split, severe in row_keys])
    cmap = LinearSegmentedColormap.from_list("gap_blue", ["#f7f9fb", "#c1d9e6"])
    left.imshow(matrix, origin="upper", aspect="auto", cmap=cmap, vmin=0, vmax=.30,
                extent=(-.5, 3.5, 5.5, -.5), interpolation="nearest")
    for row in range(6):
        for col in range(4):
            left.text(col, row, f"{matrix[row, col]:.3f}", ha="center", va="center", fontsize=8.5,
                      color="#172c39", fontweight="bold" if col == 3 else "normal")
    left.set_xticks(range(4), ["O / O", "O / P", "P / O", "P / P"])
    left.xaxis.tick_top()
    left.tick_params(axis="x", length=0, pad=8)
    left.set_yticks(range(6), [f"{split}, {severe} - 100" for split, severe in row_keys])
    left.tick_params(axis="y", length=0, pad=7)
    left.set_ylim(5.5, -.5)
    left.set_xlabel("Training / evaluation labels", labelpad=10)
    for spine in left.spines.values():
        spine.set_visible(False)
    left.set_xticks(np.arange(-.5, 4, 1), minor=True)
    left.set_yticks(np.arange(-.5, 6, 1), minor=True)
    left.grid(which="minor", color="white", linewidth=1.2)
    left.tick_params(which="minor", length=0)

    for row, (split, severe) in enumerate(row_keys):
        r = source[(source.Split == split) & (source.Severe == severe) & (source.Quantity == "paired_full_change")].iloc[0]
        color = "#176585" if severe == 11 else "#52697d"
        marker = "o" if severe == 11 else "s"
        right.errorbar(r.Gap, row, xerr=[[r.Gap - r.CI_low], [r.CI_high - r.Gap]], fmt=marker,
                       color=color, ecolor=color, elinewidth=1.2, markersize=4.5,
                       capsize=3, capthick=1, zorder=3)
        right.text(.0175, row + .29, f"{r.Gap:.3f} [{r.CI_low:.3f}, {r.CI_high:.3f}]",
                   ha="center", va="center", fontsize=7, color="#334954")
    # Short zero-reference guides leave the numeric interval labels unobstructed.
    for row in range(6):
        right.plot([0, 0], [row - .20, row + .14], color="#8b969e", linewidth=.8,
                   linestyle="--", zorder=1)
    for y in (1.5, 3.5):
        right.axhline(y, color="#e0e5e9", linewidth=.65, zorder=0)
    right.set_ylim(5.5, -.5)
    right.set_xlim(-.10, .135)
    right.set_xticks([-.10, -.05, 0, .05, .10], ["-0.10", "-0.05", "0", "0.05", "0.10"])
    right.set_yticks([])
    right.set_xlabel("Change in AUROC gap", labelpad=9)
    right.spines[["top", "right", "left"]].set_visible(False)
    right.spines["bottom"].set_color("#74808a")
    for axis, letter, title in [(left, "a", "Four-cell AUROC gap"), (right, "b", "Paired change: P/P - O/O")]:
        axis.annotate(letter, (0, 1), xycoords="axes fraction", xytext=(-11, 29), textcoords="offset points",
                      ha="left", va="bottom", fontsize=10, fontweight="bold")
        axis.annotate(title, (0, 1), xycoords="axes fraction", xytext=(0, 30), textcoords="offset points",
                      ha="left", va="bottom", fontsize=8.5, fontweight="bold")
    fig.text(.025, .067, "O: distributed labels. P: provisional source-derived labels. Panel b: paired change [95% CI].", fontsize=7.5)
    fig.text(.025, .022, "Gaps compare severe versus 100 µg/bee endpoints. Positive change means a wider gap.", fontsize=7.5)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, axes=[left, right], panel_ids=["a", "b"],
                                        json_out=QA / "Fig7_label_crossover.alignment.json",
                                        overlay_svg=QA / "Fig7_label_crossover.alignment.svg",
                                        tolerance_pt=1.5, gutter_tolerance_pt=1.5,
                                        require_panel_labels=True, strict=True)
    output = ROOT / "figures/Fig7_label_crossover"
    fig.savefig(output.with_suffix(".pdf"))
    fig.savefig(output.with_suffix(".png"), dpi=600)
    plt.close(fig)
    report = {"input": str(INPUT.relative_to(ROOT)), "input_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
              "predictions_sha256": hashlib.sha256(PRED.read_bytes()).hexdigest(), "selected_source_rows": len(source),
              "matrix_cells": 24, "paired_changes": 6, "numeric_checks": len(check_rows),
              "max_numeric_discrepancy": max(r["MaxAbsoluteDifference"] for r in check_rows),
              "inclusion_rule": "primary/keep/ECFP, all 3 splits and 2 severe endpoints, all 4 cells plus full change",
              "uncertainty": "95% scaffold-cluster percentiles, B=2000, conditional and exploratory",
              "layout": "all gap estimates on left; paired effect and uncertainty on right"}
    (QA / "FIGURE_NUMERIC_METADATA.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    caption = r"""\textbf{Matched-cohort training-by-evaluation label crossover.} \textbf{(a)} Severe-minus-broad AUROC gaps for all four combinations of distributed (O) and provisional source-derived (P) training/evaluation labels, using the primary correction rules, the keep policy and ECFP. Each row uses 828 training and 207 test compounds with identical membership across its four cells; columns show point estimates. \textbf{(b)} Change in each gap after revising both label sets (P/P minus O/O), with 95\% paired scaffold-cluster percentile intervals from 2,000 shared draws. Positive values indicate a wider gap. All three splits and both severe endpoints are shown. The intervals condition on the fitted models and implemented labels, are exploratory and are not adjusted across correction scenarios. Complete cell intervals, scores and source data accompany Supplementary Table~S23.
"""
    (QA.parent / "figure_crossover_caption.tex").write_text(caption, encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
