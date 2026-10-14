# Unseen-identity case-study protocol

This case study evaluates nine restored ECFP models on structures outside the
1,035-compound development dataset. Candidate selection was informed by
available regulatory outcomes, so the exercise is not a blinded prospective
validation.

## Models and inference

The models use the original Random, MaxMin and Time training sets and the
100, 11 and 1 microgram-per-bee labels, radius-2/1,024-bit Morgan fingerprints,
float32 Tanimoto kernels and class-balanced SVC with C=1. Model restoration
fits these nine configurations on their original training data. External
structures and outcomes do not enter fitting or model selection. The saved
restored-model scores reproduce all 1,863 corresponding benchmark test scores
to numerical precision.

Two structure-only inference batches contain eight and two candidates. They
are parts of one ten-structure case study. Structure files and predictions
are stored separately from the experimental endpoint ledger. The evaluation
joins these records by compound identity after the endpoint criteria have
been applied. All candidate structures receive predictions, including those
without eligible experimental outcomes.

## Chemical identity and endpoints

Candidate structures are screened against all development structures using
full InChIKey, connectivity, standardized parent and CAS/name evidence. A new
publication date does not establish a new or independent toxicological study.
Original source locations and study identifiers are retained where available.

The primary task is adult Apis mellifera 48-hour acute contact LD50 in
individual-dose units. Oral 48-hour observations are a separate secondary
task. Other durations remain separately identified. LD50 is distinguished
from concentration-based and chronic-dose endpoints. Inequalities and
intervals are preserved: a label is positive if LD50 is at or below the
threshold, negative if it is above the threshold, and unresolved otherwise.
No value is imputed from a prediction or an indeterminate limit observation.

## Reporting and limits

The fixed decision rule is score > 0; scores are not calibrated probabilities.
Results include each eligible endpoint, route-specific class counts, ranking
metrics, confusion matrices, paired contrasts and conditional bootstrap
intervals. Single-class ranking metrics remain undefined. Six identities
provide eleven eligible endpoints, so repeated models and endpoints do not
constitute additional independent chemicals. Small class counts, missing
severity tiers and incomplete study-identity evidence bound interpretation.
