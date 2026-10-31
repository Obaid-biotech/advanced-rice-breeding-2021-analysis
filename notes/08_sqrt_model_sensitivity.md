# Square-root response sensitivity

Run `Rscript R/08_sqrt_model_sensitivity.R` (base R only). All 159 responses were checked to be finite and nonnegative before square roots. The recorded zero is permitted and unchanged.

## Models and interpretation

Original: plant_weight ~ replication + entry. Transformed: sqrt(plant_weight) ~ replication + entry. Both predictors are factors. These remain provisional additive fixed-effects analyses of the reported 2021 RCBD, with unconfirmed units and sampling details.
A includes all 159 observations. B omits only E-46/R-3, explicitly as a sensitivity analysis, without declaring it erroneous. Original fits are refitted for comparison without replacing script 07 outputs.
Transformation changes the response scale and interpretation: the transformed model describes mean square-root grain weight, not arithmetic mean grain weight. Squaring a fitted square-root mean does not directly recover the arithmetic mean grain weight. For a transformed variable Z, E(Z^2) = [E(Z)]^2 + Var(Z); simply squaring omits the variance component. No back-transformed mean estimates are produced here.
Residual SD is in original units for the original fit and square-root units for the transformed fit; a smaller numerical SD across these scales is not evidence of a better model. Ordinary lm AIC values must not be compared directly across different response scales without accounting for the transformation of the response density. Cross-scale AIC is deliberately omitted from the exported comparison.
Model choice must not depend on obtaining a smaller p-value. The scientific response, design, sampling, and diagnostics matter. No preferred model is selected.

## Results

       scenario response_scale observations model_rank design_columns full_rank
 A_all_recorded       original          159         55             55      TRUE
 A_all_recorded           sqrt          159         55             55      TRUE
  B_omit_E46_R3       original          158         55             55      TRUE
  B_omit_E46_R3           sqrt          158         55             55      TRUE
 residual_df residual_sd_on_response_scale negative_fitted_values
         104                     7.6597548                      0
         104                     0.6529727                      0
         103                     7.6967901                      0
         103                     0.6391684                      0

 Df    Sum of Sq         RSS  F value       Pr(>F)        term       scenario
  2   966.667539  7068.53929 8.237917 4.774597e-04 replication A_all_recorded
 52 13391.639917 19493.51167 4.389355 7.266230e-11       entry A_all_recorded
  2     8.452178    52.79500 9.911710 1.148123e-04 replication A_all_recorded
 52   130.477927   174.82075 5.884962 9.620496e-15       entry A_all_recorded
  2   947.499380  7049.27888 7.997047 5.909607e-04 replication  B_omit_E46_R3
 52 12394.501713 18496.28121 4.023523 8.775621e-10       entry  B_omit_E46_R3
  2     7.208046    49.28728 8.821794 2.908361e-04 replication  B_omit_E46_R3
 52   102.217914   144.29715 4.811639 5.854578e-12       entry  B_omit_E46_R3
 response_scale
       original
       original
           sqrt
           sqrt
       original
       original
           sqrt
           sqrt

All four models are full rank. Partial F-tests use drop1(model, test = 'F'), testing each term while retaining the other. This avoids sequential term-order effects in unbalanced scenario B. The entry null is equal entry effects on the respective model's response scale after accounting for blocks; the response scales therefore have different null interpretations.

## Diagnostic layout

Each scenario's PNG has original-scale model diagnostics on the left and square-root model diagnostics on the right. Rows show residual-versus-fitted, normal Q-Q, and scale-location plots. Residuals are on each model's own response scale, not back-transformed. Only the three largest absolute standardized residuals per fit plus E-46/R-3, when present, are labelled. Smooth curves are descriptive guides.
Diagnostics cannot verify randomization or independence. Field records, layout, experimental-unit definitions, and sampling details are still needed. No assumptions are declared solved merely because a transformation was fitted.

## Visual inspection findings

Both scenario comparison PNGs were inspected. The square-root fits reduce the pronounced upward scale-location pattern seen in the original fits, especially across the central range, but some increase remains at higher fitted values. In scenario A the scale-location curve rises again at very low fitted values, where there are few observations; E-46/R-3 becomes a relatively large negative standardized residual rather than the near-zero residual seen on the original scale. Scenario B has a flatter low-to-middle fitted range after transformation, with remaining upward spread at the high end.
The square-root Q-Q plots look closer to the reference overall than the original fits, although departures and extreme residuals remain. All four fits have nonnegative fitted values in this dataset. The transformation changes which observations look unusual and has not established that constant variance or normal errors hold. Neither transformation nor omitting the zero verifies randomization, independence, or the reported experimental design.
The square-root entry F-statistics are approximately 5.885 (A) and 4.812 (B); replication F-statistics are 9.912 and 8.822. Their smaller p-values than the original fits are not a reason to select the transformation. The source zero is unchanged, and no further records are omitted.

## Outputs and preservation

outputs/sqrt_sensitivity_model_summary.csv, sqrt_sensitivity_partial_tests.csv, and sqrt_sensitivity_diagnostics.csv contain all four fits' summaries, tests, and observation diagnostics. Each scenario also has a diagnostic comparison PNG and a new models.rds containing its two fitted models.
No additional observations were deleted, lines ranked, pairwise tests run, or source observations modified. Previous scripts and outputs are preserved.

R version 4.5.2 (2025-10-31)
