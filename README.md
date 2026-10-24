# Auditing label provenance and endpoint definitions in a honey-bee toxicity benchmark

Research code, data, results and figures for the study.

| Author | Email |
|---|---|
| Yongye Huang | hiiy132321@163.com |
| Runqing Liao | 202411055222@stu.hstc.edu.cn |
| Xiaopeng Yang | 706697032@qq.com |
| Zhining Wang (corresponding author) | mathwangzhining@163.com |

All authors: School of Mathematics and Statistics, Hanshan Normal University, Chaozhou 521041, China.

The study connects three questions in the ApisTox benchmark: which tier comparisons contribute to the broad-endpoint discrimination deficit; how historical qualifier handling and aggregation determine labels; and whether source-supported partial label revisions narrow the endpoint gap. Across 15 representation-by-split settings, the broad endpoint has the lowest AUROC. Historical reconstruction identifies a qualifier-handling defect. All six primary ECFP change intervals include zero, so the tested partial revisions provide no evidence of a narrower endpoint gap.

## Download and reproduce

- [Complete research code, data, results and figures (Research.zip)](https://github.com/hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking/releases/download/research-reproducibility-2026-10-24/Research.zip)
- [Fixed reproducibility release](https://github.com/hiiiy-pm/Endpoint-Definitions-and-Label-Provenance-in-Honey-Bee-Toxicity-Benchmarking/releases/tag/research-reproducibility-2026-10-24)
- [Research overview](Research/README.md)
- [Environment and reproduction instructions](Research/docs/REPRODUCIBILITY.md)
- [Figure and table source index](Research/docs/figure_table_sources.csv)

The archive extracts into `Research/` and contains the complete computational inputs and outputs. Run the analysis commands from that directory with Python 3.12 and the pinned dependencies.

```sh
cd Research
python -m venv .venv
# Activate .venv using the command for your shell.
python -m pip install -r requirements-analysis.txt -r requirements-figures.txt
python code/run_all_analysis.py --clean
```

Run the clean analysis in a separate working copy to preserve the supplied results for comparison. External-source reconstruction, sensitivity analyses, model restoration, inference, and table and figure generation are described in the reproduction instructions.

## Research contents

| Location | Contents |
|---|---|
| `Research/code/` | Numerical analyses and the end-to-end runner |
| `Research/data/` | Benchmark records, official split membership and external-source inputs |
| `Research/results/` | Molecular representations, labels, predictions, metrics and bootstrap distributions |
| `Research/source_data/` | Machine-readable figure data |
| `Research/figures/` | Five main figures and six supplementary figures as final PDF assets |
| `Research/output/` | Source-record analyses, saved models, external evaluations and case-study outputs |
| `Research/supplementary_audits/` | Identifier screening, source evidence and sensitivity analyses |
| `Research/scripts/` | Source reconstruction and table and figure generators |
| `Research/docs/` | Data dictionary, provenance, output mapping and reproduction instructions |

## Figures

Main figures: [1](Research/figures/Figure_1.pdf), [2](Research/figures/Figure_2.pdf), [3](Research/figures/Figure_3.pdf), [4](Research/figures/Figure_4.pdf), [5](Research/figures/Figure_5.pdf).

Supplementary figures: [S1](Research/figures/Figure_S1.pdf), [S2](Research/figures/Figure_S2.pdf), [S3](Research/figures/Figure_S3.pdf), [S4](Research/figures/Figure_S4.pdf), [S5](Research/figures/Figure_S5.pdf), [S6](Research/figures/Figure_S6.pdf).

## Citation and licenses

[CITATION.cff](CITATION.cff) supplies the four-author software citation metadata. Cite the associated study and the original datasets when using these materials. Study-written code is distributed under the MIT License. Third-party data and reference code retain their original terms and attribution; see [Research/NOTICE.md](Research/NOTICE.md).
