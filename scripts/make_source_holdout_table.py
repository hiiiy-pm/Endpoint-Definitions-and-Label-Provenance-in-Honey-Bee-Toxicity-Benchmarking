from pathlib import Path
import pandas as pd
r=Path(__file__).resolve().parents[1]
m=pd.read_csv(r/'results/decomposition/08_leave_one_source_out.csv')
b=pd.read_csv(r/'results/decomposition/10_leave_one_source_out_bootstrap.csv')
lines=[r'\Needspace{12\baselineskip}\subsection*{Supplementary Table S10. Recorded-source holdout sensitivity}',r'The sources are record groups within the integrated benchmark and share its label mapping. ECFP models and the route-plus-use metadata comparator are fitted to the remaining groups. Cells report held-out AUROC; sample counts refer to compounds, not repeated predictions.',r'\begingroup\small\setlength{\tabcolsep}{4pt}',r'\begin{longtable}{llrccc}',r'\toprule Source & Model & $n$ & $\leq100$ & $\leq11$ & $\leq1$ \\ \midrule\endfirsthead',r'\toprule Source & Model & $n$ & $\leq100$ & $\leq11$ & $\leq1$ \\ \midrule\endhead',r'\bottomrule\endlastfoot']
for source in ['ECOTOX','PPDB','BPDB']:
 for model,label in [('ECFPStructure','ECFP'),('Route+Agrochemical','Route + use')]:
  d=m[(m.HeldOutSource==source)&(m.Model==model)].set_index('Threshold')
  lines.append(source+' & '+label+' & '+str(int(d.n_test.iloc[0]))+' & '+' & '.join(f'{d.loc[t,"AUROC"]:.3f}' for t in [100,11,1])+r' \\')
lines += [r'\end{longtable}',r'\Needspace{12\baselineskip}ECFP paired severe-minus-broad contrasts use scaffold-cluster simultaneous intervals for the two-contrast family.',r'\begin{longtable}{llrcc}',r'\toprule Source & Contrast & $\Delta$AUROC & 95\% lower & 95\% upper \\ \midrule\endfirsthead',r'\toprule Source & Contrast & $\Delta$AUROC & 95\% lower & 95\% upper \\ \midrule\endhead',r'\bottomrule\endlastfoot']
for source in ['ECOTOX','PPDB','BPDB']:
 for c in ['11-100','1-100']:
  d=b[(b.HeldOutSource==source)&(b.Bootstrap=='ScaffoldCluster')&(b.Comparison==c)].iloc[0]
  lines.append(source+' & '+c+' & '+f'{d.ObservedDeltaAUROC:.3f} & {d.SimultaneousCI_low:.3f} & {d.SimultaneousCI_high:.3f}'+r' \\')
lines += [r'\end{longtable}\endgroup']
(r/'output'/'latex_tables'/'supp_tables_source_holdout.tex').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print('S10 generated directly from saved source-holdout outputs; no input/result changes.')
