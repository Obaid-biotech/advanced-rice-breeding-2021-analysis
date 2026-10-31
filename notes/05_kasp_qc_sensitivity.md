# KASP QC sensitivity

Run `Rscript R/05_kasp_qc_sensitivity.R`. This script uses base R and the original audit CSVs from step 04. All input files remain unchanged.

These nine combinations are exploratory sensitivity scenarios, not validated universal thresholds. No final QC scenario is selected.

## Filtering order and denominators

1. Select eligible samples using their ORIGINAL missing-call count divided by all 502 markers. Maximum limits are inclusive (<=).
2. For each of the 502 markers, recalculate missing-call proportion using only those eligible samples. Its denominator is the retained sample count shown below, not 53.
3. Retain markers at or below the scenario's marker-missingness limit.
4. Recalculate observed allele counts from valid calls among eligible samples. Count retained markers with exactly two observed alleles as polymorphic biallelic markers.

Each scenario starts from the original audit. Sample eligibility is not recalculated after marker selection, and there is no iterative filtering or imputation. Missing calls retain the audit definition: ?, Uncallable, Bad, or blank. Unexpected strings cause the script to stop rather than silently changing that definition.
If a scenario has no eligible samples, marker missingness and allele counts are undefined (NA in the detail report), retained counts are zero, and no_retained_samples is TRUE.

## Actual comparison

| Max sample missingness | Max marker missingness | Samples | Markers | Polymorphic biallelic markers |
|---|---|---:|---:|---:|
| 10% | 10% | 42 | 452 | 417 |
| 10% | 20% | 42 | 475 | 439 |
| 10% | 30% | 42 | 487 | 449 |
| 20% | 10% | 49 | 418 | 385 |
| 20% | 20% | 49 | 467 | 433 |
| 20% | 30% | 49 | 486 | 450 |
| 30% | 10% | 52 | 414 | 381 |
| 30% | 20% | 52 | 467 | 433 |
| 30% | 30% | 52 | 484 | 448 |

Scenarios with no retained samples: 0

Excluded sample IDs (identical for all three marker limits within each sample limit):

- Sample limit 10%: NYUT_21_1;NYUT_21_4;NYUT_21_21;NYUT_21_23;NYUT_21_29;NYUT_21_31;NYUT_21_34;NYUT_21_37;NYUT_21_40;NYUT_21_41;NYUT_21_52.
- Sample limit 20%: NYUT_21_1;NYUT_21_21;NYUT_21_29;NYUT_21_37.
- Sample limit 30%: NYUT_21_37.

## Interpretation

Tighter sample limits keep fewer samples. Marker missingness can change when samples are removed, so it must be recomputed before applying the marker limit. Within a fixed sample limit, relaxing the marker limit can retain more markers.
Observed polymorphism also depends on which samples remain. If the only sample carrying an alternative allele is excluded, a formerly biallelic marker can become monomorphic. Original full-sample allele counts therefore cannot be reused. Counts here measure observed alleles, not proof that an allele is absent from the wider population.
Retained marker counts include monomorphic markers that pass missingness; polymorphic biallelic counts are a separate descriptive subset. Neither count determines the best scenario on its own.

## Files

- outputs/kasp_qc_sensitivity.csv: all nine scenarios, denominators, counts, no-sample flags, and excluded IDs.
- outputs/kasp_qc_sensitivity_marker_details.csv: recalculated missingness and allele counts for every marker in every scenario.
- outputs/kasp_qc_sensitivity.png: retained samples and polymorphic biallelic markers in separate panels. Dotted lines show original totals (53 samples and 502 markers).
No phenotype records were excluded, no marker calls were joined to phenotypes, and no marker–trait associations were tested.

R version 4.5.2 (2025-10-31)
