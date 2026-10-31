# Provisional grain weight per plant RCBD analysis

## Model and scope

Model: `plant_weight ~ replication + entry`. Both predictors are factors: categorical labels, not continuous numbers. This additive fixed-effects model accounts for replication (block) while estimating differences among these 53 entries in the reported 2021 RCBD. It does not model entry-by-block interaction separately.
plant_weight is user-confirmed grain weight per plant. Units and sampling details remain unconfirmed. recorded_yield is not analyzed here; its conversion and previously documented discrepancies remain unresolved.
Scenario A uses all 159 recorded observations including E-46/R-3 = 0. Scenario B omits only that record as a sensitivity analysis. This neither alters the source data nor declares the observation erroneous. Neither scenario is selected based on significance.

## Model checks and partial tests

       scenario observations design_columns model_rank full_rank residual_df
 A_all_recorded          159             55         55      TRUE         104
  B_omit_E46_R3          158             55         55      TRUE         103
 residual_sd
    7.659755
    7.696790

 Df  Sum of Sq       RSS      AIC  F value       Pr(>F)        term
 NA         NA  6101.872 689.9440       NA           NA      <none>
  2   966.6675  7068.539 709.3263 8.237917 4.774597e-04 replication
 52 13391.6399 19493.512 770.6203 4.389355 7.266230e-11       entry
 NA         NA  6101.779 687.2910       NA           NA      <none>
  2   947.4994  7049.279 706.0975 7.997047 5.909607e-04 replication
 52 12394.5017 18496.281 758.5113 4.023523 8.775621e-10       entry
       scenario
 A_all_recorded
 A_all_recorded
 A_all_recorded
  B_omit_E46_R3
  B_omit_E46_R3
  B_omit_E46_R3

Both fits are full rank. Residual degrees of freedom equal observations minus fitted model rank. Residual SD estimates the typical residual scale, in the unconfirmed grain-weight unit, using residual sum of squares divided by residual degrees of freedom before taking the square root.
The F-tests come from drop1(model, test = 'F'): each term is tested after accounting for the other term. This is important for the unbalanced scenario B and avoids sequential tests that depend on term order. The <none> row describes the full model, not an additional test.
The entry-effect null hypothesis is that all entry effects are equal after accounting for replication. It is a joint test, not a comparison of individual pairs. The block null hypothesis is equal replication effects after accounting for entry. Small p-values describe evidence under model assumptions; they do not establish genetic superiority or validate the design.

## Beginner explanations and diagnostic limits

lm() fits a linear model by minimizing the sum of squared residuals. A residual is observed grain weight minus its model-fitted value. Standardized residuals also account for estimated residual SD and leverage. Cook's distance summarizes how much the fitted model would change if one observation were omitted; it is not an automatic deletion rule.
Residual-versus-fitted plots look for systematic patterns; Q-Q plots compare residual quantiles with a normal reference; scale-location plots examine variation in residual spread; Cook's-distance plots highlight relative influence. Every point is labelled by entry/replication; dense labels can overlap, and CSVs provide exact values.
The usual F-tests assume independent errors with a common variance and an appropriate additive mean structure, with normal errors for exact small-sample inference. Diagnostics cannot verify randomization or independence. Field layout, randomization records, sampling details, and experimental-unit definitions must be checked separately. Plot appearances alone cannot validate those assumptions.

## Visual inspection and scenario comparison

All eight diagnostic PNGs were inspected. The patterns are similar across scenarios. Residual-versus-fitted and scale-location plots show wider residual spread at higher fitted values, raising concern about constant variance. Q-Q plots broadly follow the reference through much of the distribution but show departures, particularly in the upper tail. These are visual concerns, not formal assumption validation.
E-40/R-3 and E-51/R-2 have the largest Cook's distances in both fits, with E-51/R-1 also relatively prominent. They remain included. E-46/R-3 has a small residual and low Cook's distance in scenario A: a recorded zero need not be unusual relative to its fitted entry and block effects. This does not establish its biological validity.
Omitting E-46/R-3 reduces observations and residual df by one, slightly increases residual SD, and reduces the entry partial F-statistic from about 4.39 to 4.02. Block F changes from about 8.24 to 8.00. Both scenarios give small partial-test p-values under the stated assumptions, but neither is preferred based on those p-values. The variance pattern and unconfirmed sampling/design details qualify inference in both scenarios.

## Files and restrictions

outputs/grain_weight_rcbd_scenario_summary.csv and grain_weight_rcbd_partial_tests.csv contain the comparison. Scenario-specific diagnostic CSVs contain observed values, fitted values, residuals, standardized residuals, Cook's distances, and leverage. Each scenario has four labelled diagnostic PNGs and its own partial-test CSV.
No lines were ranked, pairwise tests run, heritability estimated, observations transformed automatically, or influential records deleted. The only omission is the explicitly requested sensitivity scenario. Existing data and scripts are preserved.

R version 4.5.2 (2025-10-31)
