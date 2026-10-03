# 最终论文研究母版

## 1. 最终研究定位

### 推荐英文题目

**Regulatory-Boundary Heterogeneity in Machine-Learning Prediction of Honey Bee Acute Toxicity: Chemical-Space Ambiguity and Out-of-Distribution Generalization**

### 推荐中文题目

**蜜蜂急性毒性机器学习预测中的监管边界异质性：化学空间模糊性与分布外泛化**

### 论文类型定位

本研究应定位为：

> **科学发现 / 机制解释型 chemoinformatics + ecotoxicology 研究**

而不是：

> 新算法/SOTA 模型论文。

核心问题不是“哪个模型 AUROC 最高”，而是：

1. 不同监管毒性边界是否具有不同的 molecular learnability？
2. 这种差异在 Random、MaxMin、Time 等不同泛化条件下是否稳定？
3. 局部 chemical neighborhood 与 scaffold organization 能否解释部分边界差异？
4. Chemical novelty 如何影响 OOD prediction reliability？

---

## 2. 数据与最终标签

原始数据：1035 个化合物，13 个字段：

`name, CID, CAS, SMILES, source, year, toxicity_type, herbicide, fungicide, insecticide, other_agrochemical, label, ppdb_level`

来源：PPDB 510、ECOTOX 441、BPDB 84。

根据已有 EPA 二分类标签与 PPDB 分级标签，无人为调阈值地重建四级风险：

| Tier | LD50 范围 | n | 解释 |
|---|---:|---:|---|
| Tier0 | >100 μg/bee | 177 | non/low toxicity |
| Tier1 | 11–100 μg/bee | 562 | moderate-low regulatory region |
| Tier2 | 1–11 μg/bee | 125 | EPA-toxic transition region |
| Tier3 | ≤1 μg/bee | 171 | highly toxic |

构建三个 cumulative endpoints：

- `Y100 = I(Tier ≥ 1)`
- `Y11 = I(Tier ≥ 2)`
- `Y1 = I(Tier ≥ 3)`

这三个任务不是互相独立的数据集，而是**同一批 molecules 上的不同 regulatory endpoints**。

---

## 3. 官方划分

三个划分全部沿用第二轮复现包中的官方重建协议：

| Split | Train | Test | 主要意义 |
|---|---:|---:|---|
| Random | 828 | 207 | IID-like |
| MaxMin | 828 | 207 | chemical-space extrapolation |
| Time | 828 | 207 | temporal extrapolation |

正式论文中不要为了 Time 结果显著而重新随机划分。

---

## 4. 分子表示

最终 robustness 采用五类表示：

1. ECFP4, 1024-bit
2. Avalon, 1024-bit
3. MACCS
4. WL-HI graph representation
5. 12 physicochemical descriptors

12 个描述符：MolWt、MolLogP、TPSA、HDonors、HAcceptors、RotBonds、RingCount、AromaticRings、FractionCSP3、HeavyAtoms、MolMR、FormalCharge。

正式 boundary statistical comparison 使用固定 ECFP/Tanimoto + class-balanced SVM，以减少“不同阈值选择不同最佳模型”造成的模型选择偏差。

---

# 5. Primary Finding 1：Regulatory-boundary heterogeneity

## 5.1 Representation robustness

5 representations × 3 splits = 15 个设置中：

**15/15 均满足：**

`AUROC(100 μg) < AUROC(11 μg)` 且 `AUROC(100 μg) < AUROC(1 μg)`。

这支持：

> 100 μg/bee boundary 的结构可学习性系统性弱于两个 severe-toxicity boundaries。

但不能写：

> 1 μg 一定严格优于 11 μg。

MaxMin/Time 中若干 representation 对 11 与 1 μg 的排序相等或反转。

## 5.2 正式 5000 次 paired bootstrap

### Molecule-level bootstrap

| Split | Comparison | ΔAUROC | 95% CI | Holm-primary p | 结论 |
|---|---|---:|---:|---:|---|
| Random | 11−100 | +0.130 | 0.009–0.256 | 0.0368 | 显著 |
| Random | 1−100 | +0.199 | 0.073–0.331 | 0.0048 | 显著 |
| MaxMin | 11−100 | +0.201 | 0.088–0.312 | 0.0004 | 显著 |
| MaxMin | 1−100 | +0.204 | 0.055–0.344 | 0.0076 | 显著 |
| Time | 11−100 | +0.061 | −0.069–0.189 | 0.3436 | 不显著 |
| Time | 1−100 | +0.101 | −0.042–0.235 | 0.3184 | 不显著 |

### Scaffold-cluster bootstrap

| Split | Comparison | ΔAUROC | 95% CI | Holm-primary p | 结论 |
|---|---|---:|---:|---:|---|
| Random | 11−100 | +0.122 | −0.018–0.243 | 0.0824 | 不显著 |
| Random | 1−100 | +0.192 | 0.054–0.313 | 0.0240 | 显著 |
| MaxMin | 11−100 | +0.199 | 0.090–0.309 | 0.0004 | 显著 |
| MaxMin | 1−100 | +0.202 | 0.060–0.337 | 0.0056 | 显著 |
| Time | 11−100 | +0.062 | −0.074–0.195 | 0.3648 | 不显著 |
| Time | 1−100 | +0.100 | −0.048–0.249 | 0.3408 | 不显著 |

### 最终表述

**可以写：**

> Boundary heterogeneity was robust across molecular representations and was statistically strongest under MaxMin chemical-space extrapolation, where severe-toxicity endpoints outperformed the 100 μg/bee endpoint by approximately 0.20 AUROC even after scaffold-cluster bootstrap.

**不能写：**

> Time split statistically confirms the same difference.

Time 仅“方向一致，但统计证据不足”。

---

# 6. Primary Finding 2：Chemical-space ambiguity

## 6.1 Neighborhood entropy

对 ECFP、Avalon、MACCS、WL-HI、Descriptors12，分别取 k=5/10/20 最近邻，计算归一化四级标签熵，并做 2000 次 label permutation。

最终 global BH-FDR：

> **15/15 representation × k permutation tests 均 q<0.05，最大 q≈0.0245。**

因此可以较强地写：

> Tier2 exhibits persistent excess local neighborhood ambiguity across molecular representations and neighborhood scales.

例如 k=10：

| Representation | Tier2 entropy | Tier2 same-tier neighbor fraction |
|---|---:|---:|
| ECFP | 0.613 | 0.213 |
| Avalon | 0.592 | 0.218 |
| MACCS | 0.627 | 0.209 |
| WL-HI | 0.627 | 0.229 |
| Descriptors12 | 0.667 | 0.166 |

注意：不能称 Tier2 “永远是 entropy 最大的唯一类别”；部分 representation 中 Tier0 或 Tier3 也可能同样高。

## 6.2 Scaffold ambiguity

主分析要求 scaffold n≥3：

| Tier | Mean purity | Mean normalized entropy | Mixed scaffold fraction |
|---|---:|---:|---:|
| Tier0 | 0.563 | 0.726 | 0.954 |
| Tier1 | 0.658 | 0.597 | 0.878 |
| Tier2 | 0.553 | 0.747 | **1.000** |
| Tier3 | 0.683 | 0.565 | 0.878 |

全局 Holm 校正：

- Tier2 vs Tier1 entropy: adjusted p = 0.0062
- Tier2 vs Tier3 entropy: adjusted p = 0.0137

在更严格 scaffold n≥5：

- Tier2 vs Tier1 仍显著：p=0.0396
- Tier2 vs Tier3 不再显著：p=0.0921

因此最终表述必须带敏感性边界：

> Tier2 scaffold ambiguity is robust relative to Tier1; the Tier2–Tier3 contrast weakens when restricting analysis to larger scaffolds.

同时 Tier0 也表现出较高 scaffold heterogeneity。因此最合理的机制叙事不是“只有 Tier2 有问题”，而是：

- Tier2：持续的 local neighborhood ambiguity / transition-region mixing；
- Tier0：较强 scaffold heterogeneity；
- 两者共同帮助解释某些监管边界的弱 structural separability。

---

# 7. Secondary Finding：Chemical novelty / applicability domain

定义：

`novelty = 1 - max_train_Tanimoto(test molecule)`

使用训练集内 5-fold OOF decision scores 做 Platt calibration，不接触 test labels；随后以 Brier loss 建立 molecule-clustered GEE：

`Brier ~ Boundary * Novelty`

## 7.1 Novelty main effect

100 μg endpoint：

- MaxMin：每 +0.1 novelty，Brier error +0.0274，95% CI 0.0124–0.0425，p=0.000365
- Time：每 +0.1 novelty，Brier error +0.0194，95% CI 0.00695–0.0319，p=0.00227

因此可以写：

> Chemical novelty was associated with deteriorating predictive reliability under distribution shift, particularly for the 100 μg/bee endpoint.

## 7.2 Boundary × Novelty interaction —— 负结果

正式 Brier-GEE 交互：0/6 contrasts 经 Holm 后显著。

最接近的是 MaxMin：

- 1 vs 100 slope difference: raw p=0.0562, Holm p=0.1124
- joint interaction p=0.144

Random joint p=0.881；Time joint p=0.452。

Binomial GEE classification-error sensitivity analysis 也不显著。

因此**绝对不能写**：

> Chemical novelty significantly moderates boundary heterogeneity.

正确写法：

> Although stratified descriptive analyses suggested larger boundary differences within familiar chemical space, formal continuous Boundary × Novelty interactions were not statistically supported.

---

# 8. Secondary analysis：Sparse physicochemical signature

使用 12 个描述符，训练集内 Elastic-Net logistic（SGD 实现）选择正则化，500 次 balanced bootstrap 计算 stability selection frequency；随后评估 Top 1/2/3/5/8/12 descriptor saturation。

主要结果：

- Random severe endpoints 可由少量 descriptors 捕获一部分信号，例如 ≤1 μg Top-3 AUROC≈0.730，Top-5≈0.757。
- MaxMin 仍有部分信号。
- Time 明显衰减；full 12 descriptors 也不足以形成稳定 temporal transfer。

因此该部分只能写为辅助发现：

> Simple physicochemical properties capture part of severe-toxicity discrimination under IID/selected OOD settings, but they do not constitute a universally transferable low-dimensional signature.

不要将 HAcceptors、LogP、TPSA 等称为“因果毒理机制”或“普适决定因子”。

---

# 9. 多重检验最终规则

必须在 Methods 中提前说明：

1. **Primary regulatory-boundary comparisons**：每个 split × bootstrap level 内对 11−100、1−100 使用 Holm correction。
2. **Neighborhood permutation family**：15 个 representation×k permutation tests 统一 BH-FDR。
3. **Neighborhood pairwise exploratory tests**：统一 BH-FDR。
4. **Scaffold predefined pairwise tests**：跨 minN=3/5 统一 Holm。
5. **Descriptor univariate associations**：36 tests 统一 BH-FDR。
6. 主要结论同时报告 effect size、95% CI、adjusted p/q，而不只报告 p-value。

---

# 10. Negative-control 与防泄漏设计

必须在论文中保留：

- 标签置换 negative control：真实模型明显优于随机标签 null；MaxMin 100 μg 最接近弱可学习边界。
- split 固定，测试集从不用于模型参数选择。
- GEE 概率校准使用训练集 OOF score 完成 Platt sigmoid。
- scaffold-cluster bootstrap 处理结构类似 analog 的聚类依赖。
- 不使用连续 LD50，因为当前冻结数据包不存在该字段；不从 categorical labels 反推。

---

# 11. 最终主结论只能保留三层

## Primary 1 — Regulatory-boundary heterogeneity

100 μg/bee endpoint 的结构可学习性稳定弱于 11 和 1 μg/bee severe endpoints，MaxMin 下证据最强。

## Primary 2 — Chemical-space ambiguity mechanism

Tier2 展现稳定的 local neighborhood ambiguity；scaffold-level mixing 提供独立支持，但 Tier0 同样具有较高 scaffold heterogeneity，说明弱可分性并非单一 Tier 所致。

## Secondary — Applicability-domain sensitivity

Chemical novelty 在 MaxMin/Time 下增加预测误差，尤其是 100 μg endpoint；但 Boundary × Novelty formal interaction 不成立。

---

# 12. 推荐论文结构

## Introduction

1. 蜜蜂急性毒性监管预测的重要性。
2. 既往研究主要关注单一阈值模型性能、QSAR、AD 或新模型。
3. 同一数据中的不同 regulatory boundaries 往往被作为等价 classification tasks，但其 structural learnability 是否相同尚不清楚。
4. 提出三个研究问题：boundary heterogeneity、chemical-space explanation、OOD reliability。

## Methods

2.1 Data source and regulatory-tier reconstruction  
2.2 Official Random/MaxMin/Time splits  
2.3 Molecular representations  
2.4 Fixed endpoint models  
2.5 Representation robustness  
2.6 Molecule/scaffold paired bootstrap  
2.7 Local neighborhood entropy and permutation tests  
2.8 Bemis–Murcko scaffold purity/entropy  
2.9 Chemical novelty and GEE reliability analysis  
2.10 Sparse physicochemical analysis  
2.11 Multiple-testing correction and negative controls

## Results

3.1 Regulatory thresholds exhibit heterogeneous learnability  
3.2 Boundary differences are strongest under MaxMin extrapolation  
3.3 Tier2 shows persistent local chemical-space ambiguity  
3.4 Scaffold organization supports toxicity-tier heterogeneity  
3.5 Chemical novelty degrades OOD reliability but does not significantly interact with boundary  
3.6 Physicochemical descriptors provide only partial low-dimensional explanation

## Discussion

1. 100 μg boundary 为什么可能缺乏清晰 structural separation。
2. Tier2 transition-region ambiguity 的含义。
3. OOD prediction 与 AD 的现实意义。
4. 为什么复杂模型不是本论文贡献重点。
5. 限制：单数据集、categorical endpoint、Time power、实验条件异质性、不能声称因果。
6. 监管风险建模应报告 boundary-specific generalizability，而不是只给一个总体 benchmark。

---

# 13. 推荐主文图表

### Figure 1 — Study design
数据 → 4 tiers → 3 regulatory endpoints → 3 splits → representations → statistical/mechanistic validation。

### Figure 2 — Representation × Split AUROC heatmap
突出 15/15 的 100 μg 弱可学习性。

### Figure 3 — Bootstrap forest plot
展示 Random/MaxMin/Time 的 ΔAUROC 和 molecule/scaffold CI。

### Figure 4 — Neighborhood ambiguity
五种 representation × k 的 Tier entropy；主图可选 k=10，k=5/20 放补充材料。

### Figure 5 — Scaffold purity/entropy
Tier0–Tier3 的 scaffold entropy/purity，n≥3 主分析，n≥5 sensitivity。

### Figure 6 — Novelty vs Brier error
展示 GEE fitted slopes；图注必须明确 interaction 不显著。

### Supplementary Figure/Table
Sparse descriptors、negative control、seen/unseen scaffold、全部 robustness matrix。

---

# 14. 投稿时最容易被攻击的点及防御

### 攻击 1：只是换了三个阈值
防御：创新不是“使用三个阈值”，而是系统研究 regulatory-boundary learnability heterogeneity，并用多 representation、OOD split、paired/scaffold bootstrap 和 chemical-space ambiguity 分析建立证据链。

### 攻击 2：AUC 差异来自模型选择
防御：正式 statistical comparison 固定同一 ECFP/Tanimoto SVM family；representation robustness 独立验证。

### 攻击 3：molecules 并不独立
防御：同时报告 molecule paired bootstrap 与 Bemis–Murcko scaffold-cluster bootstrap。

### 攻击 4：Tier2 只是样本少
防御：之前 balanced sample control 已显示 severe-vs-100 pattern 仍存在；neighborhood permutation 保持 class counts 并验证 excess mixing。

### 攻击 5：Novelty 交互是事后故事
防御：正式 continuous GEE interaction 未显著，论文明确报告负结果，不夸大。

### 攻击 6：descriptor 是机制因果
防御：descriptor analysis 明确作为 association/secondary analysis，不作 causal claim。

---

# 15. 实验冻结原则

从此版本开始：

- 不再增加新模型以追求 SOTA；
- 不再增加新的监管阈值；
- 不重新划分 Time 以追求显著；
- 不更换 novelty 定义反复寻找 interaction significance；
- 不因结果不漂亮删除 negative findings；
- 仅允许在审稿要求下进行预先说明的补充敏感性分析。

这保证后续论文写作与结果不会发生“研究者自由度”膨胀。
