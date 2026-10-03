from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]
out=root/'output'/'latex_tables'/'supp_tables.tex'
data=root/'source_data'
perf=pd.read_csv(data/'figure_primary_data.csv',dtype={'endpoint':str})
model=pd.read_csv(data/'figure_model_decomposition.csv',dtype={'endpoint':str})
matched=pd.read_csv(data/'figure_matched_n125.csv',dtype={'endpoint':str})
scaf=pd.read_csv(data/'figure_scaffold_corrected.csv')

def longtable_open(table_number, column_spec, header):
    """Repeat the table identifier and column names after each page break."""
    columns = header.count('&') + 1
    return [
        rf'\begin{{longtable}}{{{column_spec}}}',
        header + r'\endfirsthead',
        rf'\multicolumn{{{columns}}}{{l}}{{\textit{{Supplementary Table {table_number} (continued)}}}} \\',
        header + r'\endhead',
        rf'\midrule\multicolumn{{{columns}}}{{r}}{{Continued on next page}}\\\endfoot',
        r'\bottomrule\endlastfoot',
    ]


lines=[]
lines.append(r'\subsection*{Supplementary Table S1. Full test AUROC across representations, splits and endpoints}')
lines.extend(longtable_open('S1', 'llccc', r'\toprule Split & Representation & $\leq100$ & $\leq11$ & $\leq1$ \\ \midrule'))
for sp in ['Random','MaxMin','Time']:
    for rep in ['ECFP','Avalon','MACCS','WL-HI','Descriptors12']:
        d=perf[(perf['split']==sp)&(perf['representation']==rep)].set_index('endpoint')
        vals=[float(d.loc[ep,'auroc']) for ep in ['100','11','1']]
        lines.append(f'{sp} & {rep} & {vals[0]:.3f} & {vals[1]:.3f} & {vals[2]:.3f} \\\\')
    lines.append(r'\addlinespace')
lines.append(r'\end{longtable}')

lines.append(r'\subsection*{Supplementary Table S5. Test AUROC for metadata, structure and fused models}')
lines.extend(longtable_open('S5', 'llccc', r'\toprule Split & Model & $\leq100$ & $\leq11$ & $\leq1$ \\ \midrule'))
for sp in ['Random','MaxMin','Time']:
    for mdl in ['Agrochemical flags','Source + route','All metadata','ECFP structure','Structure + metadata']:
        d=model[(model['split']==sp)&(model['model']==mdl)].set_index('endpoint')
        vals=[float(d.loc[ep,'auroc']) for ep in ['100','11','1']]
        mdltex=mdl.replace('+',r'$+$')
        lines.append(f'{sp} & {mdltex} & {vals[0]:.3f} & {vals[1]:.3f} & {vals[2]:.3f} \\\\')
    lines.append(r'\addlinespace')
lines.append(r'\end{longtable}')

mb = pd.read_csv(root/'results/decomposition/06_model_decomposition_bootstrap.csv')
mb = mb[(mb.Bootstrap == 'ScaffoldCluster') & mb.PrimaryFamily]
lines.append(r'Fusion contrasts below use scaffold-cluster 95\% simultaneous intervals within each split and endpoint. The family contains structure minus metadata, fusion minus metadata and fusion minus structure; it is separate from the endpoint-contrast family. All 5,000 draws are valid for these contrasts.')
assert (mb.B_valid == 5000).all()
lines.extend(longtable_open('S5', 'llccc', r'\toprule Split & Endpoint & Contrast & $\Delta$AUROC & 95\% simultaneous CI \\ \midrule'))
model_contrasts = {'ECFPStructure-AllMetadata': 'Structure--metadata', 'ECFP+AllMetadata-AllMetadata': 'Fusion--metadata', 'ECFP+AllMetadata-ECFPStructure': 'Fusion--structure'}
for r in mb.itertuples():
    lines.append(f'{r.Split} & {r.Threshold} & {model_contrasts[r.Comparison]} & {r.ObservedDeltaAUROC:.3f} & {r.SimultaneousCI_low:.3f} to {r.SimultaneousCI_high:.3f}' + r' \\')
lines.append(r'\end{longtable}')

lines.append(r'\subsection*{Supplementary Table S9. Mean AUROC and training-subsample ranges at 125 compounds per class}')
lines.append(r'Cells give means and 2.5th--97.5th percentile ranges across 200 repeated balanced training subsamples. The test identities are fixed. These ranges describe training-subsample variability and are not independent-test confidence intervals; Figure~4a uses medians.')
lines.extend(longtable_open('S9', 'lccc', r'\toprule Split & $\leq100$ & $\leq11$ & $\leq1$ \\ \midrule'))
ms = pd.read_csv(root/'results/learning_curves/02_matched_learning_curves_summary.csv')
mc = pd.read_csv(root/'results/learning_curves/03_matched_learning_curve_contrasts.csv')
for sp in ['Random','MaxMin','Time']:
    d=ms[(ms.Split==sp)&(ms.n_per_class==125)].set_index('Threshold')
    cells=[f'{d.loc[t,"MeanAUROC"]:.3f} [{d.loc[t,"AUROC_Q025"]:.3f}, {d.loc[t,"AUROC_Q975"]:.3f}]' for t in [100,11,1]]
    lines.append(sp + ' & ' + ' & '.join(cells) + r' \\')
lines.append(r'\end{longtable}')
lines.extend(longtable_open('S9', 'llccc', r'\toprule Split & Contrast & Mean $\Delta$AUROC & 2.5th percentile & 97.5th percentile \\ \midrule'))
for r in mc[(mc.n_per_class==125)&mc.Comparison.isin(['11-100','1-100'])].itertuples():
    lines.append(f'{r.Split} & {r.Comparison} & {r.MeanDeltaAUROC:.3f} & {r.Q025:.3f} & {r.Q975:.3f}' + r' \\')
lines.append(r'\end{longtable}')

lines.append(r'\subsection*{Supplementary Table S7. Leave-one-out scaffold-entropy contrasts}')
lines.append(r'Scaffold-equal and shared-scaffold paired analyses estimate different summaries. The paired analysis restricts to scaffolds containing both compared tiers; its negative contrasts should not be interpreted as confirmation of the scaffold-equal estimand.')
sb = pd.read_csv(root/'results/geometry/08_scaffold_sensitivity_bootstrap.csv')
lines.extend(longtable_open('S7', 'llllccc', r'\toprule Min. size & Weighting & Contrast & Estimate & 95\% CI low & 95\% CI high & Scaffolds \\ \midrule'))
for _,r in scaf.iterrows():
    mn=str(r['minimum']).replace('≥',r'$\geq$')
    an='Equal' if r['analysis'] == 'Scaffold-equal' else 'Paired'
    co=str(r['contrast']).replace('−',r'$-$')
    nmin = int(str(r['minimum'])[-1])
    analysis = 'ScaffoldEqualWeighted' if r['analysis'] == 'Scaffold-equal' else 'PairedScaffoldsContainingBothTiers'
    nsc = int(sb[(sb.MinScaffoldN == nmin) & (sb.Analysis == analysis) & (sb.Comparison == str(r['contrast']).replace('−','-'))].ScaffoldsResampled.iloc[0])
    lines.append(f'{mn} & {an} & {co} & {r.estimate:.3f} & {r.low:.3f} & {r.high:.3f} & {nsc} \\\\')
lines.append(r'\end{longtable}')

boot=pd.read_csv(root/'results'/'primary'/'05_representation_paired_bootstrap.csv')
boot=boot[(boot['PrimaryComparison'])&(boot['Bootstrap']=='ScaffoldCluster')]
# Source-data exports belong to code/09_export_submission_tables.py.
# Reading this analysis table must not mutate their provenance hashes.

lines.append(r'\subsection*{Supplementary Table S2. Scaffold-cluster simultaneous intervals for the two-contrast family within each representation and split}')
assert (boot.B_valid == 5000).all()
lines.append(r'All 5,000 bootstrap draws are valid for every reported contrast (zero single-class draws). Intervals are simultaneous only within the two stated contrasts for each representation and split.')
lines.extend(longtable_open('S2', 'llcccc', r'\toprule Split & Representation & Contrast & $\Delta$AUROC & 95\% CI low & 95\% CI high \\ \midrule'))
comp_tex={'11-100':r'$\leq11-\leq100$','1-100':r'$\leq1-\leq100$'}
for sp in ['Random','MaxMin','Time']:
    for rep in ['ECFP','Avalon','MACCS','WL-HI','Descriptors12']:
        for comp in ['11-100','1-100']:
            d=boot[(boot['Split']==sp)&(boot['Representation']==rep)&(boot['Comparison']==comp)].iloc[0]
            reptex=rep.replace('WL-HI','WL--HI')
            lines.append(
                f'{sp} & {reptex} & {comp_tex[comp]} & {d.ObservedDeltaAUROC:.3f} & '
                f'{d.SimultaneousCI_low:.3f} & {d.SimultaneousCI_high:.3f} \\\\'
            )
    lines.append(r'\addlinespace')
lines.append(r'\end{longtable}')

rev = root / 'results' / 'sensitivity'
EPL = {100: r'$\leq100$', 11: r'$\leq11$', 1: r'$\leq1$'}

fam = pd.read_csv(rev / '08d_ten_contrast_family_intervals.csv')
fam = fam[fam.Bootstrap == 'ScaffoldCluster']
lines.append(r'\subsection*{Supplementary Table S3. Scaffold-cluster simultaneous intervals over the ten contrasts within each split}')
lines.extend(longtable_open('S3', 'llccc', r'\toprule Split & Representation & Contrast & $\Delta$AUROC & 95\% CI \\ \midrule'))
for sp in ['Random', 'MaxMin', 'Time']:
    for _, r in fam[fam.Split == sp].iterrows():
        c = comp_tex[r.Comparison]
        lines.append(f'{sp} & {r.Representation.replace("WL-HI", "WL--HI")} & {c} & '
                     f'{r.Delta:.3f} & {r.CI_low:.3f} to {r.CI_high:.3f} \\\\')
    lines.append(r'\addlinespace')
lines.append(r'\end{longtable}')

lrn = pd.read_csv(rev / '08e_learner_sensitivity.csv')
lines.append(r'\subsection*{Supplementary Table S4. Test AUROC for additional learners on ECFP4 bits}')
lines.append(r'Point-estimate sensitivity. Logistic regression selects $C$ by five-fold training-only AUROC; the selected values at 100/11/1 are 0.10/0.10/0.10 for Random and 0.01/0.10/0.10 for MaxMin and Time. Random forest uses 500 trees and \texttt{balanced\_subsample} weights; details are in Supplementary Table~S21a.')
lines.extend(longtable_open('S4', 'llccc', r'\toprule Split & Learner & $\leq100$ & $\leq11$ & $\leq1$ \\ \midrule'))
names = {'LogisticL2': 'L2 logistic regression', 'RandomForest': 'Random forest'}
for sp in ['Random', 'MaxMin', 'Time']:
    for ln in ['LogisticL2', 'RandomForest']:
        d = lrn[(lrn.Split == sp) & (lrn.Learner == ln)].set_index('Endpoint').AUROC
        lines.append(f'{sp} & {names[ln]} & {d[100]:.3f} & {d[11]:.3f} & {d[1]:.3f} \\\\')
    lines.append(r'\addlinespace')
lines.append(r'\end{longtable}')

cx = pd.read_csv(rev / '08a_cross_endpoint_evaluation.csv')
lines.append(r'\subsection*{Supplementary Table S12. Cross-endpoint evaluation (rows: endpoint the model was trained for; columns: endpoint whose labels were used for evaluation)}')
lines.append(r'Cells give AUROC and 95\% scaffold-cluster percentile intervals from 2,000 draws. Each off-diagonal entry reuses the row model\textquotesingle s score vector; no model is refitted for the column labels.')
lines.extend(longtable_open('S12', 'lllccc', r'\toprule Split & Representation & Trained for & $\leq100$ & $\leq11$ & $\leq1$ \\ \midrule'))
for sp in ['Random', 'MaxMin', 'Time']:
    for rep in ['ECFP', 'WL-HI']:
        for tr_ep in [100, 11, 1]:
            d = cx[(cx.Split == sp) & (cx.Representation == rep) & (cx.TrainEndpoint == tr_ep)]
            d = d.set_index('EvalEndpoint')
            cells = [f'{d.AUROC[e]:.3f} ({d.CI_low[e]:.2f}--{d.CI_high[e]:.2f})' for e in (100, 11, 1)]
            lines.append(f'{sp} & {rep.replace("WL-HI", "WL--HI")} & {EPL[tr_ep]} & ' + ' & '.join(cells) + r' \\')
        lines.append(r'\addlinespace')
lines.append(r'\end{longtable}')

mm = pd.read_csv(rev / '08c_matched_preprocessing_metadata.csv')
lines.append(r'\subsection*{Supplementary Table S6. Metadata-only, structure-only and fused models under identical standardization}')
lines.extend(longtable_open('S6', 'llcccc', r'\toprule Split & Endpoint & Metadata & Metadata without route & ECFP & Fused (all / without route) \\ \midrule'))
for sp in ['Random', 'MaxMin', 'Time']:
    for ep in [100, 11, 1]:
        a = mm[(mm.Split == sp) & (mm.Endpoint == ep) & (mm.MetadataSet == 'AllMetadata')].iloc[0]
        b = mm[(mm.Split == sp) & (mm.Endpoint == ep) & (mm.MetadataSet == 'MetadataNoRoute')].iloc[0]
        lines.append(f'{sp} & {EPL[ep]} & {a.MetadataAUROC:.3f} & {b.MetadataAUROC:.3f} & '
                     f'{a.StructureAUROC:.3f} & {a.FusedAUROC:.3f} / {b.FusedAUROC:.3f} \\\\')
    lines.append(r'\addlinespace')
lines.append(r'\end{longtable}')

bs = pd.read_csv(rev / '08b_brier_skill.csv')
lines.append(r'\subsection*{Supplementary Table S13. Brier score of calibrated ECFP forecasts against a constant training-prevalence forecast}')
lines.extend(longtable_open('S13', 'llcccc', r'\toprule Split & Endpoint & Training prevalence & Model Brier & Constant Brier & Brier skill \\ \midrule'))
for sp in ['Random', 'MaxMin', 'Time']:
    for ep in [100, 11, 1]:
        r = bs[(bs.Split == sp) & (bs.Endpoint == ep)].iloc[0]
        lines.append(f'{sp} & {EPL[ep]} & {r.TrainingPrevalence:.3f} & {r.ModelBrier:.4f} & '
                     f'{r.ConstantBrier:.4f} & {r.BrierSkill:.3f} \\\\')
    lines.append(r'\addlinespace')
lines.append(r'\end{longtable}')

out.write_text('\n'.join(lines), encoding='utf-8', newline='\n')
