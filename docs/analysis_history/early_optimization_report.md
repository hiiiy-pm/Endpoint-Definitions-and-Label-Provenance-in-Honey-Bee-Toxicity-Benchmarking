# ApisTox 第二轮优化结果

## 1. 优化原则

本轮只使用训练集内部完成模型选择，重点验证：

1. MaxMin 外层使用 repeated MaxMin inner validation；
2. Time 外层使用 forward temporal inner validation；
3. 引入图结构信息，而不是继续增加相似 fingerprint；
4. 对 positive class weight 和 decision threshold 做训练内选择；
5. 保留第一轮 Equal-MK-WKELM 作为固定基线；
6. 评价 AUROC、AUPRC、MCC、Sensitivity、Specificity，而不是只追求单一指标。

## 2. 第一轮基线

| Split   |   AUROC |   AUPRC |    MCC |   Sensitivity |   Specificity |
|:--------|--------:|--------:|-------:|--------------:|--------------:|
| Random  |  0.8697 |  0.778  | 0.6674 |        0.678  |        0.9459 |
| MaxMin  |  0.831  |  0.6202 | 0.4279 |        0.4524 |        0.9273 |
| Time    |  0.7645 |  0.5786 | 0.47   |        0.439  |        0.9518 |

## 3. 第二轮模型

| Split   | Model                           |   AUROC |   AUPRC |    MCC |   Sensitivity |   Specificity |
|:--------|:--------------------------------|--------:|--------:|-------:|--------------:|--------------:|
| MaxMin  | WL-HI SVM (shift-aware)         |  0.8344 |  0.5781 | 0.4527 |        0.5952 |        0.8727 |
| MaxMin  | Current-MK + WL-HI ensemble     |  0.843  |  0.5971 | 0.4344 |        0.5952 |        0.8606 |
| MaxMin  | Avalon SVM (high-MCC candidate) |  0.7723 |  0.5788 | 0.4823 |        0.4524 |        0.9515 |
| Time    | WL-HI SVM (temporal-aware)      |  0.7936 |  0.5492 | 0.4579 |        0.3902 |        0.9639 |
| Time    | ECFP+WL WKELM (temporal-aware)  |  0.7771 |  0.6237 | 0.4754 |        0.3902 |        0.9699 |

## 4. 关键发现

### MaxMin

第一轮基线：
- AUROC = 0.8310
- MCC = 0.4279
- Sensitivity = 0.4524

WL-HI SVM：
- AUROC = 0.8344
- MCC = 0.4527
- Sensitivity = 0.5952

Current-MK + WL-HI ensemble：
- AUROC = 0.8430
- MCC = 0.4344
- Sensitivity = 0.5952

Avalon SVM：
- AUROC = 0.7723
- MCC = 0.4823

这说明 MaxMin 当前存在明显的 Pareto trade-off：
- 图结构模型/ensemble 更有利于 AUROC 和 toxic sensitivity；
- Avalon SVM 更有利于 MCC；
- 暂时没有单一模型同时支配所有指标。

### Time

第一轮基线：
- AUROC = 0.7645
- MCC = 0.4700

WL-HI SVM：
- AUROC = 0.7936
- MCC = 0.4579

ECFP+WL WKELM：
- AUROC = 0.7771
- AUPRC = 0.6237
- MCC = 0.4754

Time split 的 AUROC 已得到明显改善，说明图结构信息对于时间外推具有价值。

## 5. 学术解释

本轮最重要的结果不是“所有指标同时大幅上涨”，而是确认：
- Random-CV 最优模型并不等价于 shift-robust 模型；
- Time split 对 temporal-aware validation 有响应；
- MaxMin 需要更强的结构互补信息；
- AUROC、MCC、Sensitivity 之间存在真实的多目标权衡；
- 后续应采用 Pareto/多目标学习，而不是为了某一个 test 指标继续事后调参。

## 6. 重要限制

WL-HI 是本研究实现的 Weisfeiler-Lehman histogram-intersection / optimal-assignment-style graph kernel，
并非声称与 2026 benchmark 的 WL-OA 实现逐行一致。
正式投稿前应复现公开 benchmark 的原始 WL-OA 代码作为强基线。

另外，本轮属于模型开发阶段。由于研究过程中已经查看过官方 test 结果，
正式论文应在最终方法锁定后补充严格的 repeated nested validation，
并避免继续根据 test 结果反复改模型。
