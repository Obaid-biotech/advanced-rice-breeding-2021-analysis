# Candidate evidence summary

Run `Rscript R/12_candidate_evidence_summary.R` (base R only). Existing scripts and output structures were inspected. Prior data and outputs are preserved. Validated key joins use the user-confirmed phenotype E-n to NYUT_21_n correspondence; this is not independent field-record verification.

## E-42 and E-45: all 502 original markers

Valid calls in both: 487
          category count
    matching_valid   487
 mismatching_valid     0
  missing_E42_only     2
  missing_E45_only     1
      missing_both    12
Matching and mismatching valid counts partition valid-in-both; do not add valid-in-both a second time. The five disjoint categories total exactly 502. Missing calls remain ?, Bad, Uncallable, or blank; canonical reversed pairs agree, while heterozygous and homozygous genotypes differ. All mismatches and their original calls are exported, even if the resulting CSV contains headers only.

    kasp_id sample_name master_plate master_well phenotype_entry
 NYUT_21_42     RC21075     TSK14060         B06            E-42
 NYUT_21_45     RC21080     TSK14060         E06            E-45
Laboratory metadata are extracted read-only from Sheet1 A:D via the workbook ZIP/XML using base R and checked against the mapping. Different sample names or wells document different labels/locations, not independently established biological identity.

Original phenotype values side by side (plant_weight means grain weight per plant, unit unconfirmed):

 replication E-42_plant_height E-42_tillers E-42_panicle_length
         R-1          96.66667     12.33333            25.37778
         R-2          97.66667     24.00000            25.87778
         R-3          97.66667     16.66667            24.88889
 E-42_filled_grains E-42_unfilled_grains E-42_grain_weight_100
          106.33333             9.222222              2.370000
          117.88889            18.066667              2.466667
           91.66667            12.444444              2.603333
 E-42_plant_weight E-42_recorded_yield E-45_plant_height E-45_tillers
          42.77333            5810.757          95.33333     21.66667
          61.26000            8322.171         100.00000     16.33333
          27.72333            3766.215          97.00000     14.00000
 E-45_panicle_length E-45_filled_grains E-45_unfilled_grains
            28.36667           127.8889             14.11111
            25.76667           123.0000             15.22222
            26.01111           135.5556             12.11111
 E-45_grain_weight_100 E-45_plant_weight E-45_recorded_yield
              2.360000          60.14333            8170.472
              2.403333          53.73000            7299.221
              2.593333          34.24000            4651.504
recorded_yield is a separate unresolved measurement. The 135.85 relationship is unverified; audit 02 discrepancies are carried into the evidence table and are not corrected.

## Possible explanations and needed evidence

Possible explanations include distinct materials not distinguished by this targeted marker set, shared ancestry or closely related material, a seed/sample handling or labelling issue, or a genotyping/data transfer issue. These are alternatives to investigate, not conclusions. Neither identical observed panel calls nor different phenotype measurements establishes biological identity or its absence.
Review field entry books and randomization maps, seed inventories and pedigree/accession records, collection labels and chain of custody, DNA extraction logs, plate maps, assay cluster/raw fluorescence files, and export/transcription history. Independently collect traceable plants or seed lots, repeat DNA extraction and genotyping with controls, and consider additional informative genome-distributed markers or sequencing. Replicate sampling within entries may help distinguish within-entry variation from handling problems. No identities are corrected and no crosses are recommended.

## Ten-candidate evidence

- E-3: original A grain-weight estimate 48.127 (pointwise 95% CI 39.357 to 56.896); missingness 6.97%, heterozygosity among valid calls 17.13%. Supports exploratory follow-up within this trial only. Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.
- E-10: original A grain-weight estimate 46.033 (pointwise 95% CI 37.264 to 54.803); missingness 4.18%, heterozygosity among valid calls 0.62%. Supports exploratory follow-up within this trial only. Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.
- E-11: original A grain-weight estimate 43.374 (pointwise 95% CI 34.605 to 52.144); missingness 3.39%, heterozygosity among valid calls 0.41%. Supports exploratory follow-up within this trial only. Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.
- E-15: original A grain-weight estimate 47.227 (pointwise 95% CI 38.457 to 55.996); missingness 7.17%, heterozygosity among valid calls 1.93%. Supports exploratory follow-up within this trial only. Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.
- E-16: original A grain-weight estimate 45.396 (pointwise 95% CI 36.626 to 54.165); missingness 3.78%, heterozygosity among valid calls 1.04%. Supports exploratory follow-up within this trial only. Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.
- E-37: original A grain-weight estimate 46.088 (pointwise 95% CI 37.318 to 54.857); missingness 30.88%, heterozygosity among valid calls 17.00%. Supports exploratory follow-up within this trial only. Specific cautions: exceeds exploratory 10% sample missingness limit (retained). Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.
- E-40: original A grain-weight estimate 48.876 (pointwise 95% CI 40.106 to 57.645); missingness 15.94%, heterozygosity among valid calls 0.71%. Supports exploratory follow-up within this trial only. Specific cautions: recorded_yield factor discrepancy requires source review; exceeds exploratory 10% sample missingness limit (retained). Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.
- E-42: original A grain-weight estimate 43.919 (pointwise 95% CI 35.149 to 52.689); missingness 2.79%, heterozygosity among valid calls 0.41%. Supports exploratory follow-up within this trial only. Specific cautions: Combined versus replication recorded_yield discrepancy requires source review; E-42/E-45 observed-call agreement requires identity provenance review, without an identity conclusion. Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.
- E-45: original A grain-weight estimate 49.371 (pointwise 95% CI 40.601 to 58.141); missingness 2.59%, heterozygosity among valid calls 0.41%. Supports exploratory follow-up within this trial only. Specific cautions: E-42/E-45 observed-call agreement requires identity provenance review, without an identity conclusion. Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.
- E-51: original A grain-weight estimate 58.137 (pointwise 95% CI 49.367 to 66.906); missingness 2.59%, heterozygosity among valid calls 1.64%. Supports exploratory follow-up within this trial only. Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.

candidate_evidence_summary.csv contains membership in each model's five-entry set; original-scale estimates and provisional mean-response pointwise CIs for A and B; square-root estimates and CIs on their own scale; descriptive means of the other six traits; a separately labelled unresolved recorded_yield mean; original molecular QC metrics and exploratory flags; and specific audit counts.
Other trait means are descriptions, not favorable directions. No direction, weight, or composite score is assigned without a defined breeding objective. All ten remain retained, including E-37. Molecular missingness uses all 502 markers; heterozygosity uses valid calls. Selection on these same trial data may exaggerate performance. Agreement across the four related fits is sensitivity evidence, not independent validation or probability of superiority.
A includes all observations, B omits only E-46/R-3 as previously defined sensitivity. Square-root means are not squared into arithmetic means. Confidence intervals are pointwise, model-dependent, not simultaneous, not prediction intervals, and not adjusted for candidate selection. Residual variance concerns persist. Diagnostics cannot verify randomization or independence; units and sampling details remain unresolved. These results concern this trial, not broad adaptation, and further trials are needed.

## Outputs

outputs/candidate_evidence_summary.csv; e42_e45_marker_category_counts.csv; e42_e45_marker_total_validation.csv; e42_e45_all_marker_comparison.csv; e42_e45_mismatching_markers.csv; e42_e45_laboratory_metadata.csv; e42_e45_phenotypes_side_by_side.csv.
No candidates removed, identities corrected, marker associations tested, significance letters assigned, or crosses recommended.

R version 4.5.2 (2025-10-31)
