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
    for mdl in ['Agrochemical flags','All metadata','ECFP structure','Structure + metadata']:
        d=model[(model['split']==sp)&(model['model']==mdl)].set_index('endpoint')
        vals=[float(d.loc[ep,'auroc']) for ep in ['100','11','1']]
        mdltex=mdl.replace('+',r'$+$')
        lines.append(f'{sp} & {mdltex} & {vals[0]:.3f} & {vals[1]:.3f} & {vals[2]:.3f} \\\\')
    lines.append(r'\addlinespace')
lines.append(r'\end{longtable}')

lines.append(r'\subsection*{Supplementary Table S9. Matched training sets at 125 compounds per class}')
lines.append(r'\begin{tabularx}{\linewidth}{lCCC}')
lines.append(r'\toprule Split & $\leq100$ & $\leq11$ & $\leq1$ \\ \midrule')
for sp in ['Random','MaxMin','Time']:
    d=matched[matched['split']==sp].set_index('endpoint')
    vals=[float(d.loc[ep,'auroc']) for ep in ['100','11','1']]
    lines.append(f'{sp} & {vals[0]:.3f} & {vals[1]:.3f} & {vals[2]:.3f} \\\\')
lines.append(r'\bottomrule\end{tabularx}')

lines.append(r'\subsection*{Supplementary Table S7. Leave-one-out scaffold-entropy contrasts}')
lines.extend(longtable_open('S7', 'llllcc', r'\toprule Minimum size & Weighting & Contrast & Estimate & 95\% CI low & 95\% CI high \\ \midrule'))
for _,r in scaf.iterrows():
    mn=str(r['minimum']).replace('≥',r'$\geq$')
    an=str(r['analysis']).replace('-',r'--')
    co=str(r['contrast']).replace('−',r'$-$')
    lines.append(f'{mn} & {an} & {co} & {r.estimate:.3f} & {r.low:.3f} & {r.high:.3f} \\\\')
lines.append(r'\end{longtable}')

boot=pd.read_csv(root/'results'/'primary'/'05_representation_paired_bootstrap.csv')
boot=boot[(boot['PrimaryComparison'])&(boot['Bootstrap']=='ScaffoldCluster')]
# Source-data exports belong to code/09_export_submission_tables.py.
# Reading this analysis table must not mutate their provenance hashes.

lines.append(r'\subsection*{Supplementary Table S2. Scaffold-cluster paired endpoint contrasts for all representations}')
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
