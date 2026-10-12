# Label-crossover figure QA

The final PDF and 600 dpi PNG were inspected after the final layout change. The figure is 7.16 by 3.8 inches. Both panels retain the same six row centers and plot-area top/bottom edges. Render-time alignment passed at 1.5 pt tolerance. The minimum rendered PDF font size is 7 pt; PDF text remains editable. The final collision audit found no failures or warnings.

| Panel | Evidence | Estimate and uncertainty | Unit | Visual result |
|---|---|---|---|---|
| a | All four training/evaluation cells for three splits and two severe endpoints | Twenty-four gap point estimates; marginal cell intervals remain in the source CSV, while panel b shows the paired uncertainty that addresses the figure question | 207 test compounds per split, with identical members across the four cells | Every cell is labelled; subdued blue encodes gap magnitude, and bold P/P values locate the full-revision cell. The 24 text overlays are intentionally contained in heatmap cells. No overlap, clipping or omitted result was observed. |
| b | Full-revision gap minus original gap | Six paired changes, with 95% scaffold-cluster percentile intervals from 2,000 draws | Same test-compound/scaffold draws across label cells and endpoints; fitted models fixed | Numeric changes and intervals are printed below each point. Circle/square markers also identify the two endpoint contrasts. All intervals include zero; no significance symbols are used. Short zero guides leave the numeric labels unobstructed. |

The whole-figure read confirms that panel a supplies the four-cell context and panel b carries the paired inference. The caption defines the conditional, exploratory uncertainty and sample counts. A darker heatmap value means a larger gap, not better performance. The plot retains all requested cells and all six changes.

Static-preflight warnings were reviewed: the requested submission figure is vector PDF with a PNG preview, so a TIFF is unnecessary; the random-number generator only reproduces bootstrap weights for a numeric verification check and does not create simulated observations. The static tool's SVG-export inference refers to the QA alignment overlay, not an additional scientific SVG export.

Numeric verification recomputed every gap from saved compound scores and regenerated all six paired-change intervals using the original seeds and draw order. Thirty checks passed, with maximum discrepancy below 1e-16. The source table records 30 entries (24 cells and six paired changes) and includes the source marginal intervals and membership hashes.
