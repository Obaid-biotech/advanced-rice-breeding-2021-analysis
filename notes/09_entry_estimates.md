# Entry estimates with equal block weighting

Run `Rscript R/09_entry_estimates.R` (base R). Model definitions, factor levels, and data checks follow inspected scripts 07 and 08: plant_weight ~ replication + entry and sqrt(plant_weight) ~ replication + entry. A includes all 159 observations; B omits only E-46/R-3 as sensitivity, without declaring the zero erroneous. All four models are full rank.

## Estimand and calculation

For each entry, predict at R-1, R-2, and R-3 and average with weights 1/3 each. Construct those three design-matrix rows using the fitted model's terms, factor levels, and contrasts, then average them to obtain L. Estimate = L beta; SE = sqrt(L vcov(model) L'). The full covariance matrix, including off-diagonal covariances, is used. Averaging three individual prediction SEs would be incorrect.
Equal weighting targets the same three-block average for every entry, even E-46 in B. Its R-3 response is model-estimated, not an imputed observation inserted into the dataset. An adjusted mean is a model-based average response, not an individual plant observation. Units and sampling details remain unconfirmed.
Pointwise 95% confidence limits are estimate +/- qt(0.975, residual df) * SE. These are mean-response confidence intervals, not prediction intervals. Pointwise means each interval has its stated coverage under the assumptions when considered separately; the collection does not have simultaneous 95% coverage. Overlapping or nonoverlapping intervals are not a formal pairwise test.
Square-root estimates and intervals stay on the modeled square-root scale. They are not squared or labelled arithmetic mean grain weight. recorded_yield and its unresolved conversion are not analyzed.

## Verification

All four models have 53 estimates in E-1 to E-53 numeric order. Estimates are checked against averaged predict() values. SEs are checked independently using the design matrix QR decomposition. Scenario A original-scale adjusted means agree with raw entry means within 1e-10, as expected for the balanced additive model.
       scenario response_scale observations rank residual_df
 A_all_recorded       original          159   55         104
 A_all_recorded           sqrt          159   55         104
  B_omit_E46_R3       original          158   55         103
  B_omit_E46_R3           sqrt          158   55         103
 max_prediction_check_error max_se_qr_check_error
               7.105427e-15          3.552714e-15
               8.881784e-16          3.330669e-16
               7.105427e-15          2.664535e-15
               0.000000e+00          2.775558e-16
 max_raw_mean_check_error_A_original
                        1.385558e-13
                                  NA
                                  NA
                                  NA

## Sensitivity changes

The following table identifies the largest absolute A-to-B changes, not a breeding or performance ranking. Differences concern means on each respective response scale and are not tested for significance.
E-46 is the only entry whose estimated mean changes beyond numerical tolerance (1e-10). Its original-scale increase is about 0.1252 and its square-root-scale increase about 0.6201. All other changes are floating-point noise (less than 2e-13), not meaningful differences between scenarios; their apparent order in the five-row table has no substantive interpretation. Balanced entries retain their three-block raw means under the additive model. Interval widths can still change when model variance and degrees of freedom change.
 response_scale entry estimate_A estimate_B change_B_minus_A absolute_change
       original  E-46  3.5833333  3.7085208    1.2518745e-01   1.2518745e-01
       original   E-7 28.9622222 28.9622222   -1.5631940e-13   1.5631940e-13
       original   E-1 36.9222222 36.9222222   -9.9475983e-14   9.9475983e-14
       original  E-44 32.4677778 32.4677778   -7.1054274e-14   7.1054274e-14
       original   E-6 34.0488889 34.0488889   -6.3948846e-14   6.3948846e-14
           sqrt  E-46  1.5292093  2.1493061    6.2009682e-01   6.2009682e-01
           sqrt   E-7  5.3217520  5.3217520   -1.4210855e-14   1.4210855e-14
           sqrt   E-9  4.2581737  4.2581737    6.2172489e-15   6.2172489e-15
           sqrt  E-13  5.1605114  5.1605114   -6.2172489e-15   6.2172489e-15
           sqrt   E-6  5.8337339  5.8337339   -5.3290705e-15   5.3290705e-15

E-46 estimates and intervals:
       scenario response_scale entry  estimate standard_error residual_df
 A_all_recorded       original  E-46 3.5833333     4.42236150         104
 A_all_recorded           sqrt  E-46 1.5292093     0.37699396         104
  B_omit_E46_R3       original  E-46 3.7085208     5.45986836         103
  B_omit_E46_R3           sqrt  E-46 2.1493061     0.45340661         103
    lower_95   upper_95
 -5.18637548 12.3530421
  0.78161612  2.2768025
 -7.11984036 14.5368819
  1.25008104  3.0485312

## Interpretation and limits

Intervals rely on the additive model, independent errors, common residual variance, and normal-error assumptions for exact t inference. Residual variance concerns from scripts 07 and 08 persist, so these are provisional model-based intervals, not assumption-free uncertainty estimates. Diagnostics cannot verify randomization or independence.
The original-scale E-46 intervals extend below zero because the unconstrained linear-model t interval is symmetric. The limits are not clipped or interpreted as physically possible negative grain weights. Both scenario-comparison plots were inspected; A and B estimates largely coincide except for E-46, whose interval is wider in B. The square-root plot remains entirely on the modeled scale.
These estimates concern the reported 2021 trial and these blocks, not broad adaptation across environments or years. They do not establish genetic superiority. Neither scenario nor response scale is selected based on significance. No pairwise tests, significance-group letters, or breeding shortlist is produced.

## Outputs

outputs/entry_estimates_all_models.csv contains all 212 estimates; four scenario/scale CSVs each contain 53 rows. entry_estimates_averaged_design_rows.csv records the L rows; entry_estimates_validation.csv records numerical checks. entry_estimates_scenario_changes.csv retains all changes in numeric order within each scale, and entry_estimates_largest_sensitivity_changes.csv reports the five largest absolute changes per scale. Two scenario_comparison PNGs show A and B estimates and pointwise intervals in numeric entry order.
All existing scripts, model files, and source observations are preserved.

R version 4.5.2 (2025-10-31)
