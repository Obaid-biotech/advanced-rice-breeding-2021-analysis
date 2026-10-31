# Run from the project folder: Rscript R/07_grain_weight_rcbd.R
# Base R provisional fixed-effects analysis of the reported 2021 RCBD.
# plant_weight is user-confirmed grain weight per plant; units and sampling
# details remain unconfirmed. Source values are never modified.
phenotype <- read.csv("outputs/phenotype_2021_imported.csv", stringsAsFactors = FALSE)
stopifnot(all(c("entry", "replication", "plant_weight") %in% names(phenotype)),
          is.numeric(phenotype$plant_weight), nrow(phenotype) == 159,
          !anyNA(phenotype[c("entry", "replication", "plant_weight")]),
          all(is.finite(phenotype$plant_weight)),
          !anyDuplicated(phenotype[c("entry", "replication")]))
stopifnot(setequal(phenotype$entry, paste0("E-", 1:53)),
          setequal(phenotype$replication, c("R-1", "R-2", "R-3")),
          all(table(phenotype$entry, phenotype$replication) == 1))
# Factors tell R these labels are categories, not numeric trends.
phenotype$entry <- factor(phenotype$entry, levels = paste0("E-", 1:53))
phenotype$replication <- factor(phenotype$replication, levels = c("R-1", "R-2", "R-3"))
omit <- phenotype$entry == "E-46" & phenotype$replication == "R-3"
stopifnot(sum(omit) == 1, phenotype$plant_weight[omit] == 0)
cat("Model: plant_weight ~ replication + entry\n",
    "Additive fixed effects for blocks and these 53 entries.\n",
    "Partial entry null: all entry effects are equal after accounting for blocks.\n",
    "Units and sampling details are unconfirmed.\n\n", sep = "")
dir.create("outputs", showWarnings = FALSE)
dir.create("notes", showWarnings = FALSE)
scenarios <- list(A_all_recorded = phenotype, B_omit_E46_R3 = phenotype[!omit, ])
summaries <- list()
tests <- list()

for (scenario in names(scenarios)) {
  data <- scenarios[[scenario]]
  # lm() fits the additive model by least squares. No interaction is fitted:
  # there is one observation per entry/block cell in the full dataset.
  model <- lm(plant_weight ~ replication + entry, data = data, na.action = na.fail)
  design <- model.matrix(model)
  full_rank <- model$rank == ncol(design)
  if (!full_rank || df.residual(model) <= 0) stop("Model is rank deficient or lacks residual df.")
  summaries[[scenario]] <- data.frame(
    scenario = scenario, observations = nobs(model), design_columns = ncol(design),
    model_rank = model$rank, full_rank = full_rank,
    residual_df = df.residual(model), residual_sd = sigma(model)
  )
  # drop1() compares the full model with each term omitted while retaining
  # the other term. This avoids order-dependent sequential tests in scenario B.
  partial <- as.data.frame(drop1(model, test = "F"))
  partial$term <- rownames(partial)
  rownames(partial) <- NULL
  partial$scenario <- scenario
  tests[[scenario]] <- partial
  write.csv(partial, paste0("outputs/grain_weight_rcbd_", scenario, "_partial_tests.csv"), row.names = FALSE)
  diagnostics <- data.frame(
    scenario = scenario, entry = data$entry, replication = data$replication,
    observed = data$plant_weight, fitted = fitted(model), residual = residuals(model),
    standardized_residual = rstandard(model), cooks_distance = cooks.distance(model),
    leverage = hatvalues(model)
  )
  write.csv(diagnostics, paste0("outputs/grain_weight_rcbd_", scenario, "_diagnostics.csv"), row.names = FALSE)

  # Every plotted observation is labelled by entry/replication. Large PNGs
  # and small labels allow readers to zoom in; overlapping labels are possible.
  labels <- paste(data$entry, data$replication, sep = "/")
  red <- data$entry == "E-46" & data$replication == "R-3"
  colors <- ifelse(red, "#BC3427", "#305F73")
  start_plot <- function(kind, title) {
    png(paste0("outputs/grain_weight_rcbd_", scenario, "_", kind, ".png"),
        width = 2400, height = 1800, res = 180)
    par(mar = c(5, 5, 5, 2))
    title
  }
  label_points <- function(x, y, cols = colors, labs = labels) {
    text(x, y, labels = labs, pos = 3, offset = 0.35, cex = 0.48, col = cols)
  }
  title <- start_plot("residual_vs_fitted", paste(scenario, "— residuals versus fitted"))
  plot(diagnostics$fitted, diagnostics$residual, pch = 16, col = colors,
       xlab = "Fitted grain weight per plant (unit unconfirmed)",
       ylab = "Residual (observed - fitted)", main = title)
  abline(h = 0, lty = 2, col = "gray50")
  lines(lowess(diagnostics$fitted, diagnostics$residual), col = "#B5771B", lwd = 2)
  label_points(diagnostics$fitted, diagnostics$residual)
  dev.off()
  title <- start_plot("normal_qq", paste(scenario, "— normal Q-Q"))
  ord <- order(diagnostics$standardized_residual)
  theoretical <- qnorm(ppoints(nrow(data)))
  sorted <- diagnostics$standardized_residual[ord]
  plot(theoretical, sorted, pch = 16, col = colors[ord],
       xlab = "Theoretical normal quantiles", ylab = "Standardized residuals", main = title)
  qqline(diagnostics$standardized_residual, col = "#B5771B", lwd = 2)
  label_points(theoretical, sorted, colors[ord], labels[ord])
  dev.off()
  title <- start_plot("scale_location", paste(scenario, "— scale-location"))
  spread <- sqrt(abs(diagnostics$standardized_residual))
  plot(diagnostics$fitted, spread, pch = 16, col = colors,
       xlab = "Fitted grain weight per plant (unit unconfirmed)",
       ylab = expression(sqrt(abs("Standardized residual"))), main = title)
  lines(lowess(diagnostics$fitted, spread), col = "#B5771B", lwd = 2)
  label_points(diagnostics$fitted, spread)
  dev.off()
  title <- start_plot("cooks_distance", paste(scenario, "— Cook's distance"))
  plot(seq_len(nrow(data)), diagnostics$cooks_distance, type = "h", col = colors,
       xlab = "Observation index in scenario (labels identify entry/block)",
       ylab = "Cook's distance", main = title)
  points(seq_len(nrow(data)), diagnostics$cooks_distance, pch = 16, col = colors)
  label_points(seq_len(nrow(data)), diagnostics$cooks_distance)
  dev.off()
}
summary_table <- do.call(rbind, summaries)
rownames(summary_table) <- NULL
partial_tests <- do.call(rbind, tests)
rownames(partial_tests) <- NULL
write.csv(summary_table, "outputs/grain_weight_rcbd_scenario_summary.csv", row.names = FALSE)
write.csv(partial_tests, "outputs/grain_weight_rcbd_partial_tests.csv", row.names = FALSE)
print(summary_table, row.names = FALSE)
print(partial_tests, row.names = FALSE)

# Notes contain reproducible numerical results. Diagnostic interpretation is
# recorded separately after visual inspection, without changing either fit.
notes <- c(
  "# Provisional grain weight per plant RCBD analysis", "",
  "## Model and scope", "",
  "Model: `plant_weight ~ replication + entry`. Both predictors are factors: categorical labels, not continuous numbers. This additive fixed-effects model accounts for replication (block) while estimating differences among these 53 entries in the reported 2021 RCBD. It does not model entry-by-block interaction separately.",
  "plant_weight is user-confirmed grain weight per plant. Units and sampling details remain unconfirmed. recorded_yield is not analyzed here; its conversion and previously documented discrepancies remain unresolved.",
  "Scenario A uses all 159 recorded observations including E-46/R-3 = 0. Scenario B omits only that record as a sensitivity analysis. This neither alters the source data nor declares the observation erroneous. Neither scenario is selected based on significance.", "",
  "## Model checks and partial tests", "",
  capture.output(print(summary_table, row.names = FALSE)), "",
  capture.output(print(partial_tests, row.names = FALSE)), "",
  "Both fits are full rank. Residual degrees of freedom equal observations minus fitted model rank. Residual SD estimates the typical residual scale, in the unconfirmed grain-weight unit, using residual sum of squares divided by residual degrees of freedom before taking the square root.",
  "The F-tests come from drop1(model, test = 'F'): each term is tested after accounting for the other term. This is important for the unbalanced scenario B and avoids sequential tests that depend on term order. The <none> row describes the full model, not an additional test.",
  "The entry-effect null hypothesis is that all entry effects are equal after accounting for replication. It is a joint test, not a comparison of individual pairs. The block null hypothesis is equal replication effects after accounting for entry. Small p-values describe evidence under model assumptions; they do not establish genetic superiority or validate the design.", "",
  "## Beginner explanations and diagnostic limits", "",
  "lm() fits a linear model by minimizing the sum of squared residuals. A residual is observed grain weight minus its model-fitted value. Standardized residuals also account for estimated residual SD and leverage. Cook's distance summarizes how much the fitted model would change if one observation were omitted; it is not an automatic deletion rule.",
  "Residual-versus-fitted plots look for systematic patterns; Q-Q plots compare residual quantiles with a normal reference; scale-location plots examine variation in residual spread; Cook's-distance plots highlight relative influence. Every point is labelled by entry/replication; dense labels can overlap, and CSVs provide exact values.",
  "The usual F-tests assume independent errors with a common variance and an appropriate additive mean structure, with normal errors for exact small-sample inference. Diagnostics cannot verify randomization or independence. Field layout, randomization records, sampling details, and experimental-unit definitions must be checked separately. Plot appearances alone cannot validate those assumptions.", "",
  "## Visual inspection and scenario comparison", "",
  "All eight diagnostic PNGs were inspected. The patterns are similar across scenarios. Residual-versus-fitted and scale-location plots show wider residual spread at higher fitted values, raising concern about constant variance. Q-Q plots broadly follow the reference through much of the distribution but show departures, particularly in the upper tail. These are visual concerns, not formal assumption validation.",
  "E-40/R-3 and E-51/R-2 have the largest Cook's distances in both fits, with E-51/R-1 also relatively prominent. They remain included. E-46/R-3 has a small residual and low Cook's distance in scenario A: a recorded zero need not be unusual relative to its fitted entry and block effects. This does not establish its biological validity.",
  "Omitting E-46/R-3 reduces observations and residual df by one, slightly increases residual SD, and reduces the entry partial F-statistic from about 4.39 to 4.02. Block F changes from about 8.24 to 8.00. Both scenarios give small partial-test p-values under the stated assumptions, but neither is preferred based on those p-values. The variance pattern and unconfirmed sampling/design details qualify inference in both scenarios.", "",
  "## Files and restrictions", "",
  "outputs/grain_weight_rcbd_scenario_summary.csv and grain_weight_rcbd_partial_tests.csv contain the comparison. Scenario-specific diagnostic CSVs contain observed values, fitted values, residuals, standardized residuals, Cook's distances, and leverage. Each scenario has four labelled diagnostic PNGs and its own partial-test CSV.",
  "No lines were ranked, pairwise tests run, heritability estimated, observations transformed automatically, or influential records deleted. The only omission is the explicitly requested sensitivity scenario. Existing data and scripts are preserved.", "",
  R.version.string
)
writeLines(notes, "notes/07_grain_weight_rcbd.md")
