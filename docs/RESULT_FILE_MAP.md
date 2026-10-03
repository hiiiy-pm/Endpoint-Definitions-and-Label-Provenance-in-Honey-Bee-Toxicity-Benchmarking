# Result file map

| Result | Source file |
|---|---|
| Dataset, Tier and split validation | `results/final/00_data_validation.json`; `results/decomposition/00_split_endpoint_counts.csv` |
| Tier composition table/heatmap | `results/decomposition/01_tier_composition_long.csv`; `02_tier_composition_wide.csv`; `03_tier_composition_association_tests.csv` |
| Five-representation AUROC/AP (legacy column `AUPRC`) | `results/primary/01_representation_performance.csv` |
| Primary endpoint scores | `results/primary/02_test_predictions.csv` |
| Molecule/scaffold paired endpoint bootstrap | `results/primary/03_ecfp_paired_bootstrap.csv`; distributions in `.npz` |
| Same bootstrap for all five representations | `results/primary/05_representation_paired_bootstrap.csv`; distributions in `.npz` |
| Cross-representation primary-contrast support counts | `results/primary/06_cross_representation_primary_support.csv` |
| Endpoint-ordering count | `results/primary/04_endpoint_ordering_robustness.csv` |
| Metadata/structure/combined performance | `results/decomposition/04_model_decomposition_performance.csv` |
| Decomposition predictions | `results/decomposition/05_model_decomposition_predictions.csv` |
| Decomposition paired bootstrap | `results/decomposition/06_model_decomposition_bootstrap.csv`; distributions in `.npz` |
| Insecticide-stratified performance | `results/decomposition/07_insecticide_stratified_performance.csv` |
| Leave-one-source-out sensitivity | `results/decomposition/08_leave_one_source_out.csv`; `09_leave_one_source_out_predictions.csv`; `10_leave_one_source_out_bootstrap.csv` |
| Matched learning curves | `results/learning_curves/01_matched_learning_curves_long.csv`; `02_matched_learning_curves_summary.csv`; `03_matched_learning_curve_contrasts.csv` |
| Neighborhood entropy summary | `results/geometry/01_neighborhood_entropy_summary.csv` |
| Global and conditional permutations | `results/geometry/02_conditional_neighborhood_permutation.csv`; null distributions in `.npz`; stratum sizes in `03_permutation_stratum_sizes.csv` |
| Corrected scaffold descriptors/context | `results/geometry/04_scaffold_level_descriptors.csv`; `05_leave_one_out_scaffold_context.csv`; `06_scaffold_equal_weight_tier_context.csv`; `07_scaffold_context_summary.csv` |
| Scaffold sensitivity bootstrap | `results/geometry/08_scaffold_sensitivity_bootstrap.csv`; distributions in `.npz` |
| Calibrated predictions and novelty | `results/applicability/01_calibrated_predictions_and_novelty.csv` |
| Boundary×Novelty contrasts | `results/applicability/02_boundary_novelty_interactions.csv` |
| Boundary-specific Brier slopes | `results/applicability/03_boundary_specific_brier_slopes.csv` |
| Joint interaction tests | `results/applicability/04_joint_interaction_tests.csv` |
| Sparse descriptor supplement | `results/supplementary/01_elasticnet_stability_selection.csv`; `02_topk_descriptor_performance.csv`; `03_descriptor_consensus.csv` |
| Main figures | `figures/Fig1`–`Fig6` PDF/PNG/SVG/TIFF; supplementary `figures/FigS1`–`FigS4` |
| Final key numbers and claim matrix | `results/final/01_key_numbers.json`; `09_final_evidence_matrix_english.csv`; `FINAL_RESULTS_SUMMARY_ENGLISH.md` |
| Descriptor sensitivity (Supplementary Table S8) and complete fitting audit | `results/supplementary/02_topk_descriptor_performance.csv`; `output/descriptor_audit/` |
| Source concordance (Supplementary Table S15) | `output/external_data_audit/label_concordance_summary.csv`; row-level concordance and source observations |
| Main Table 5 / Supplementary Tables S18a-b | `output/external_evaluation/frozen_metrics_ECFP.csv`; `frozen_contrasts_ECFP.csv`; `frozen_compound_members.csv` |
| Supplementary Fig. S4 | `source_data/figure_external_reassessment.csv`; `scripts/make_label_source_figure.py` |
| Unseen identities (Supplementary Tables S19a-b) | `output/unseen_identity_case_study/data_curation/frozen_scoreable_48h_endpoints.csv`; complete endpoint ledger |
| Case-study results (Supplementary Table S20) | `output/unseen_identity_case_study/evaluation/pilot_performance.csv`; case scores and contrasts |
| Independent numerical checks | `verification/verify_round1.py`; `verification/verify_round2.py` |
| Final export chain | `code/09_export_submission_tables.py`; `results/final/EXPORT_PROVENANCE.json` |
| Tier-pair AUROC (Fig. 4b; Supplementary Table S11) | `results/tier_boundary/01_tier_pair_auroc.csv`; `source_data/figure_tier_pair_auroc.csv` |
| Tier1 exclusion (Fig. 4c; Table 3; Supplementary Table S14) | `results/tier_boundary/02_tier1_exclusion_auroc.csv`; `02b_tier1_refit_test_predictions.csv`; `03_tier1_exclusion_decomposition.csv`; `source_data/figure_tier1_exclusion.csv` |
| Concordance by tier and record source (Table 4; Fig. 5a; Supplementary Table S16) | `results/tier_boundary/04_concordance_by_tier.csv`; `05_concordance_by_record_source.csv`; `source_data/figure_concordance_by_tier.csv` |
| Original consensus qualifier tracing (Supplementary Table S17) | `results/tier_boundary/06_ecotox_qualifier_propagation.csv`; `07_ecotox_qualifier_by_tier.csv`; `QUALIFIER_PROPAGATION_SUMMARY.json` |
| Factorial rule comparison (Fig. 5b; Supplementary Table S22) | `results/label_audit/03_rule_transition_counts.csv`; `source_data/figure_label_rule_comparison.csv` |
| Key numbers of the tier-boundary analyses | `results/tier_boundary/TIER_BOUNDARY_SUMMARY.json` |
| Main-text Tables 3-5 rows | `scripts/make_main_tables.py` (writes `output/main_tables/`) |
| Figure QA | `output/figure_qa/` (alignment and collision reports for Fig1-Fig6 and FigS1-FigS4) |
| Supplementary tables (LaTeX) | `output/latex_tables/supplementary_tables_ordered.tex`, assembled by `scripts/assemble_supplementary_tables.py` |
| Complete rule labels / heterogeneity / upstream cut-off selection | `results/label_audit/02_compound_labels_by_rule.csv`; `04_determining_group_heterogeneity.csv`; `05_upstream_cutoff_rules.csv` |
| Correction crossover (main Table 6; Supplementary Table S23) | `results/label_audit/06_crossover_auroc.csv`; `07_crossover_gap.csv`; `08_crossover_label_changes.csv`; `09_crossover_predictions.csv.gz`; `CROSSOVER_SUMMARY.json` |
| Non-nested extreme stress test | `results/label_audit/LABEL_NESTING_CHECK.csv`; `LABEL_NESTING_CHECK.json` |
| Record cases, source flow, overlap and metadata limits (Supplementary Table S24) | `results/record_audit/`; `docs/RECORD_CASE_AUDIT.md` |
| Independent revision and full-pipeline verification | `results/reproducibility/`; `docs/CACHE_VALIDATION.md`; `docs/END_TO_END_RERUN.md` |
| Analysis-history evidence | `docs/ANALYSIS_HISTORY.md`; `docs/analysis_history/EVIDENCE_MANIFEST.json` |
