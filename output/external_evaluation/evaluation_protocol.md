# Saved-score reassessment protocol

This retrospective analysis asks whether AUROC(11) minus AUROC(100) and
AUROC(1) minus AUROC(100) persist when the same held-out scores are evaluated
against independently retrieved source records. External-source labels and
aggregate disagreement counts were inspected during analysis development;
this is not a prospective preregistration.

## Sources and cohorts

OpenFoodTox contact candidate records and EPA S1 adult-contact records with
S4 source mappings supply the regulatory observations. Eligibility uses
chemical identity, the official split membership and saved held-out scores.
Identity checks include full InChIKey, connectivity, standardized parent,
CAS and name conflicts. Labels are aggregated within structure and route by
unanimous determinate interval labels. Paired endpoint comparisons use one
common cohort with all three labels determined.

The broad OpenFoodTox screen, its 48-hour/material-field screen and the EPA
cohort are reported separately. Original ApisTox contact-route subsets provide
the more direct source-label comparison. Different original routes also imply
endpoint-route mismatch. Same-route agreement does not establish matching
duration, formulation or underlying study. The EPA sensitivity screen requires
no recorded source SMILES/InChI disagreement and a nonmissing MRID; it is a
conservative screen rather than a diagnosis of every excluded identity.

## Models and statistics

ECFP is primary; four other saved representations provide descriptive
robustness analyses. Each Random, MaxMin and Time model is evaluated on its
own eligible set. There is no fitting, tuning, calibration or outcome-based
threshold selection in this reassessment. Overlapping test compounds are not
new chemical test cases, and predictions across splits or sources are not
treated as independent observations.

Original and regulatory-source labels are compared on exactly the same
compounds and scores. The full 207-compound test performance is descriptive
context, not a different-cohort comparator for attributing relabeling effects.
Outputs include AUROC, average precision, class counts and undefined metrics
for single-class cases. ECFP uncertainty uses 5,000 compound bootstrap draws
paired across labels and thresholds, with valid and degenerate counts.
The percentile intervals are exploratory. Cases with fewer than five
compounds in a constituent class are marked sparse; a numerical interval
does not resolve that limited information.

## Interpretation

Retention, attenuation and reversal of the endpoint gap are all reported.
Source disagreement cannot be assigned exclusively to censoring or label
error without matched original-study evidence. Database separation does not
establish independent experimental studies. New candidate identities are
curated and evaluated separately; predicted values never supply ground truth.
