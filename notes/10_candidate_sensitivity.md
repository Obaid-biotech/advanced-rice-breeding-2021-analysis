# Exploratory candidate sensitivity

Run `Rscript R/10_candidate_sensitivity.R` (base R). This is exploratory prioritization by grain weight per plant, not a declaration of statistically or genetically superior lines.
The requested set sizes 5 and 10 are explicit exploratory choices, not biologically validated cutoffs. Within each of four models, entries are ordered by descending mean estimate. The kth estimate is the cutoff; include every estimate >= cutoff - 1e-10. This absolute tolerance is on that model's response scale and addresses numerical ties, not biological equivalence. Numeric entry ID only breaks display ties, never membership ties. Actual set sizes are reported.

## Candidate sets

- A_all_recorded_original_5 (actual n=5): E-51, E-45, E-40, E-3, E-15
- A_all_recorded_original_10 (actual n=10): E-51, E-45, E-40, E-3, E-15, E-37, E-10, E-16, E-42, E-11
- A_all_recorded_sqrt_5 (actual n=5): E-51, E-45, E-40, E-3, E-15
- A_all_recorded_sqrt_10 (actual n=10): E-51, E-45, E-40, E-3, E-15, E-37, E-16, E-10, E-42, E-11
- B_omit_E46_R3_original_5 (actual n=5): E-51, E-45, E-40, E-3, E-15
- B_omit_E46_R3_original_10 (actual n=10): E-51, E-45, E-40, E-3, E-15, E-37, E-10, E-16, E-42, E-11
- B_omit_E46_R3_sqrt_5 (actual n=5): E-51, E-45, E-40, E-3, E-15
- B_omit_E46_R3_sqrt_10 (actual n=10): E-51, E-45, E-40, E-3, E-15, E-37, E-16, E-10, E-42, E-11

## Overlaps

 requested_size       scenario original_size sqrt_size overlap_count
              5 A_all_recorded             5         5             5
              5  B_omit_E46_R3             5         5             5
             10 A_all_recorded            10        10            10
             10  B_omit_E46_R3            10        10            10
 union_count                                      overlap_ids original_only
           5                          E-51;E-45;E-40;E-3;E-15              
           5                          E-51;E-45;E-40;E-3;E-15              
          10 E-51;E-45;E-40;E-3;E-15;E-37;E-10;E-16;E-42;E-11              
          10 E-51;E-45;E-40;E-3;E-15;E-37;E-10;E-16;E-42;E-11              
 sqrt_only
          
          
          
          

Requested size 5: 5 entries in all four models: E-3, E-15, E-40, E-45, E-51
Requested size 10: 10 entries in all four models: E-3, E-10, E-11, E-15, E-16, E-37, E-40, E-42, E-45, E-51

## Estimates, intervals, and molecular QC

The complete comparison CSV contains all 53 entries for each requested size, four membership indicators and actual sizes, and all four means and provisional pointwise intervals. Candidate details are also saved separately. Original-scale means target arithmetic mean grain weight equally averaged over blocks; square-root estimates target mean square-root grain weight. They are different quantities. Square-root estimates and limits are not squared or labelled arithmetic means. Units and sampling details remain unconfirmed.
Candidate phenotype IDs are joined by explicit keys through sample_id_mapping.csv to kasp_sample_audit.csv. Unique IDs, full ID coverage, and agreement of sample names are checked. The phenotype-to-KASP correspondence remains user-confirmed, not independently field-verified. Every candidate's KASP missingness is missing calls divided by all 502 markers; heterozygosity is heterozygous calls divided by valid calls. The QC metrics come from the original unfiltered KASP audit. No candidate is excluded on molecular QC grounds.
A uses all recorded observations; B omits only E-46/R-3 as sensitivity without declaring the observation erroneous. Agreement across models is sensitivity evidence, not independent validation and not a probability of superiority: all fits share the same trial and nearly all observations.
Selecting candidates from the same data can exaggerate apparent performance (selection optimism). Intervals are pointwise and model-dependent, not simultaneous or selection-adjusted. Residual variance concerns remain; interval overlap is not a formal pairwise test. These results concern this trial, not broad adaptation, and further trials are needed.

## Outputs

- outputs/candidate_sensitivity_comparison.csv: complete wide comparison in numeric ID order within each requested size.
- outputs/candidate_sensitivity_candidate_details.csv: union of candidates, with estimates, intervals, mapping, and molecular QC.
- outputs/candidate_model_membership.csv: complete long membership report with cutoffs and tolerance.
- outputs/candidate_ordered_estimates.csv: descending estimates within each model.
- outputs/candidate_scale_overlaps.csv: original versus square-root overlap within each scenario.
- outputs/candidate_sensitivity_membership.png: filled membership circles, in numeric entry order.
No source data or previous results were changed. No marker associations, significance letters, or genetic-superiority labels are produced. No final set or model is selected.

R version 4.5.2 (2025-10-31)
