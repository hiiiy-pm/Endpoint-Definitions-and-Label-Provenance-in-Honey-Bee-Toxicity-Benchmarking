"""Build the identifier-screen SI insert directly from the saved sensitivity tables."""
from __future__ import annotations
import csv
from pathlib import Path

OUT = Path(__file__).resolve().parent


def read(name):
    with (OUT / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    auc = read("current_identifier_flagged_exclusion_auroc.csv")
    gaps = read("current_identifier_flagged_exclusion_contrasts.csv")
    lookup = {(r["Split"], r["Stage"], int(r["Threshold"])): r for r in auc}
    gap_lookup = {(r["Split"], r["Stage"], int(r["Severe"])): r for r in gaps}
    stages = [("original_all", "Original/all"), ("frozen_reduced", "Frozen/reduced"), ("refit_reduced", "Refit/reduced")]
    parts = [r"""The screening set was fixed before the new model fits and contained 33 current identifier-consistency flags, comprising six fragment-parent-related differences and 27 remaining differences requiring adjudication. These are not 33 confirmed historical errors. We removed these identities from both partitions of each official split, retaining the distributed labels and frozen molecular representations for every remaining compound. No SMILES, year or chemical-use flag was replaced. The retained training/test counts were 804/198 (Random), 811/191 (MaxMin) and 808/194 (Time).

The ECFP/Tanimoto SVM used $C=1$ and balanced class weights recalculated within each retained training set. Nine original-model fits first reproduced saved scores with maximum absolute difference $2.22\times10^{-16}$. We then evaluated three stages: the original saved scores on all test identities (Original/all); the same saved scores on retained test identities (Frozen/reduced); and nine refitted models on those retained test identities (Refit/reduced). Thus the evaluation-only and refitting changes remain separately inspectable.

For all stages and endpoints, 5,000 shared draws resampled scaffold clusters on each original test cohort; reduced stages retained only their eligible observations within each draw. Intervals for the two severe-minus-broad contrasts are simultaneous centred-bootstrap intervals within a split and stage. They condition on the fitted models, fixed screening set and labels. All 5,000 draws were jointly valid in each family. Intervals on original points were regenerated under this shared sensitivity resampling and do not replace the primary-analysis intervals.

\begin{center}
\footnotesize
\begin{tabularx}{\linewidth}{llcXXX}
\toprule
Split & Stage & Train/test $n$ & AUROC, 100 & AUROC, 11 & AUROC, 1 \\
\midrule"""]
    for split in ("Random", "MaxMin", "Time"):
        for stage, label in stages:
            rows = [lookup[split, stage, t] for t in (100, 11, 1)]
            assert len({(r["n_train"], r["n_test"]) for r in rows}) == 1
            scores = [f'{float(r["AUROC"]):.3f} ({r["n_pos"]}/{r["n_neg"]})' for r in rows]
            parts.append(" & ".join([split, label, f'{rows[0]["n_train"]}/{rows[0]["n_test"]}', *scores]) + r" \\")
    parts.append(r"""\bottomrule
\end{tabularx}
\end{center}
Numbers in parentheses are positive/negative test counts at the indicated cut-off in $\mu$g per bee. Original and Frozen refer to the same saved score vectors, with different evaluation membership.

\begin{center}
\footnotesize
\begin{tabularx}{\linewidth}{llXX}
\toprule
Split & Stage & 11 minus 100 (simultaneous 95\% CI) & 1 minus 100 (simultaneous 95\% CI) \\
\midrule""")
    for split in ("Random", "MaxMin", "Time"):
        for stage, label in stages:
            rows = [gap_lookup[split, stage, t] for t in (11, 1)]
            values = [f'{float(r["Gap"]):.3f} ({float(r["SimultaneousCI_low"]):.3f} to {float(r["SimultaneousCI_high"]):.3f})' for r in rows]
            parts.append(" & ".join([split, label, *values]) + r" \\")
    parts.append(r"""\bottomrule
\end{tabularx}
\end{center}
The broad endpoint remained lowest in all three refitted ECFP evaluations. The MaxMin refit contrasts were 0.212 (simultaneous 95\% CI 0.062--0.362) and 0.209 (0.060--0.359). All paired intervals for total changes in the endpoint gaps included zero. These selected-cohort deletion results do not adjudicate individual structures or establish prospective performance. Complete per-compound scores, class counts, pointwise AUROC intervals and paired evaluation/refit components are in the \path{supplementary_audits/statistics/current_identifier_flagged_exclusion_*} files.""")
    (OUT / "si_identifier_screen_insert.tex").write_text("\n".join(parts) + "\n", encoding="utf-8")
    print(OUT / "si_identifier_screen_insert.tex")


if __name__ == "__main__":
    main()
