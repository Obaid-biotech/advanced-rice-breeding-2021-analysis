# Run from the project folder: Rscript R/09_entry_estimates.R
# Base R. Model definitions and scenarios follow validated scripts 07 and 08.
d <- read.csv("outputs/phenotype_2021_imported.csv", stringsAsFactors = FALSE)
entries <- paste0("E-", 1:53)
blocks <- c("R-1", "R-2", "R-3")
stopifnot(all(c("entry", "replication", "plant_weight") %in% names(d)),
          nrow(d) == 159, is.numeric(d$plant_weight),
          !anyNA(d[c("entry", "replication", "plant_weight")]),
          all(is.finite(d$plant_weight)), all(d$plant_weight >= 0),
          setequal(d$entry, entries), setequal(d$replication, blocks),
          all(table(d$entry, d$replication) == 1))
d$entry <- factor(d$entry, levels = entries)
d$replication <- factor(d$replication, levels = blocks)
omit <- d$entry == "E-46" & d$replication == "R-3"
stopifnot(sum(omit) == 1, d$plant_weight[omit] == 0)
scenarios <- list(A_all_recorded = d, B_omit_E46_R3 = d[!omit, ])
dir.create("outputs", showWarnings = FALSE)
dir.create("notes", showWarnings = FALSE)
estimates <- list()
checks <- list()
design_rows <- list()
for (scenario in names(scenarios)) {
  models <- list(
    original = lm(plant_weight ~ replication + entry, data = scenarios[[scenario]], na.action = na.fail),
    sqrt = lm(sqrt(plant_weight) ~ replication + entry, data = scenarios[[scenario]], na.action = na.fail)
  )
  for (scale in names(models)) {
    model <- models[[scale]]
    key <- paste(scenario, scale, sep = "_")
    stopifnot(model$rank == ncol(model.matrix(model)), df.residual(model) > 0)
    # Three prediction rows per entry, one for each block. Equal weights
    # include R-3 for E-46 in B even though that observed cell is omitted.
    grid <- expand.grid(replication = blocks, entry = entries, stringsAsFactors = FALSE)
    grid$entry <- factor(grid$entry, levels = entries)
    grid$replication <- factor(grid$replication, levels = blocks)
    X <- model.matrix(delete.response(terms(model)), grid,
                      contrasts.arg = model$contrasts, xlev = model$xlevels)
    X <- X[, names(coef(model)), drop = FALSE]
    # L has one averaged design row per entry, in numeric entry order.
    L <- t(vapply(entries, function(e) colMeans(X[grid$entry == e, , drop = FALSE]),
                  numeric(ncol(X))))
    colnames(L) <- colnames(X)
    beta <- coef(model)
    V <- vcov(model)[names(beta), names(beta), drop = FALSE]
    estimate <- as.vector(L %*% beta)
    # Var(L beta) = L V L'; off-diagonal covariances are essential.
    variances <- diag(L %*% V %*% t(L))
    stopifnot(all(variances >= 0))
    se <- sqrt(variances)
    critical <- qt(0.975, df = df.residual(model))
    result <- data.frame(
      scenario = scenario, response_scale = scale, entry = entries,
      estimate = estimate, standard_error = se, residual_df = df.residual(model),
      lower_95 = estimate - critical * se, upper_95 = estimate + critical * se,
      block_weight_R1 = 1/3, block_weight_R2 = 1/3, block_weight_R3 = 1/3
    )
    # Independent numerical checks: averaging predictions agrees with L beta;
    # fitting through the design QR agrees with the full covariance calculation.
    prediction <- predict(model, newdata = grid)
    prediction_means <- vapply(entries, function(e) mean(prediction[grid$entry == e]), 0)
    QR <- qr(model.matrix(model))
    R <- qr.R(QR)
    Lp <- L[, QR$pivot, drop = FALSE]
    via_qr <- sigma(model) * sqrt(rowSums(t(backsolve(R, t(Lp), transpose = TRUE))^2))
    raw_difference <- NA_real_
    if (scenario == "A_all_recorded" && scale == "original") {
      raw_means <- vapply(entries, function(e) mean(d$plant_weight[d$entry == e]), 0)
      raw_difference <- max(abs(estimate - raw_means))
      stopifnot(raw_difference < 1e-10)
    }
    stopifnot(max(abs(estimate - prediction_means)) < 1e-10,
              max(abs(se - via_qr)) < 1e-10, nrow(result) == 53)
    checks[[key]] <- data.frame(scenario = scenario, response_scale = scale,
      observations = nobs(model), rank = model$rank, residual_df = df.residual(model),
      max_prediction_check_error = max(abs(estimate - prediction_means)),
      max_se_qr_check_error = max(abs(se - via_qr)),
      max_raw_mean_check_error_A_original = raw_difference)
    estimates[[key]] <- result
    design_rows[[key]] <- data.frame(scenario = scenario, response_scale = scale,
                                    entry = entries, L, check.names = FALSE)
    write.csv(result, paste0("outputs/entry_estimates_", key, ".csv"), row.names = FALSE)
  }
}
all_estimates <- do.call(rbind, estimates)
validation <- do.call(rbind, checks)
write.csv(all_estimates, "outputs/entry_estimates_all_models.csv", row.names = FALSE)
write.csv(validation, "outputs/entry_estimates_validation.csv", row.names = FALSE)
write.csv(do.call(rbind, design_rows), "outputs/entry_estimates_averaged_design_rows.csv", row.names = FALSE)
changes <- list()
largest <- list()
for (scale in c("original", "sqrt")) {
  a <- estimates[[paste0("A_all_recorded_", scale)]]
  b <- estimates[[paste0("B_omit_E46_R3_", scale)]]
  stopifnot(identical(a$entry, b$entry))
  change <- data.frame(response_scale = scale, entry = entries,
    estimate_A = a$estimate, estimate_B = b$estimate,
    change_B_minus_A = b$estimate - a$estimate,
    absolute_change = abs(b$estimate - a$estimate))
  changes[[scale]] <- change
  # Largest sensitivity changes, not ranking entries by performance.
  largest[[scale]] <- change[head(order(change$absolute_change, decreasing = TRUE), 5), ]
  # Three panels preserve numeric order and give each interval enough room.
  png(paste0("outputs/entry_estimates_", scale, "_scenario_comparison.png"),
      width = 2500, height = 1900, res = 180)
  par(mfrow = c(1, 3), mar = c(5, 4, 3, 1), oma = c(2, 0, 4, 0))
  xr <- range(c(a$lower_95, a$upper_95, b$lower_95, b$upper_95))
  xr <- xr + c(-0.03, 0.03) * diff(xr)
  for (indices in list(1:18, 19:36, 37:53)) {
    y <- rev(seq_along(indices))
    plot(NA, xlim = xr, ylim = c(0.4, length(indices) + 0.6), yaxt = "n",
         xlab = if (scale == "original") "Mean grain weight (unit unconfirmed)" else "Mean sqrt(grain weight) (unit unconfirmed)",
         ylab = "", main = paste(entries[min(indices)], "to", entries[max(indices)]))
    abline(h = y, col = "gray93")
    axis(2, at = y, labels = entries[indices], las = 1, cex.axis = 0.85)
    for (which_scenario in c("A", "B")) {
      z <- if (which_scenario == "A") a else b
      offset <- if (which_scenario == "A") 0.14 else -0.14
      color <- if (which_scenario == "A") "#286C87" else "#C05A25"
      segments(z$lower_95[indices], y + offset, z$upper_95[indices], y + offset, col = color, lwd = 1.6)
      points(z$estimate[indices], y + offset, pch = if (which_scenario == "A") 16 else 17,
             col = color, cex = 0.8)
    }
  }
  mtext(paste(if (scale == "original") "Original response" else "Square-root response",
              "— equal-block entry means and pointwise 95% CIs"), side = 3, outer = TRUE, line = 2, cex = 1.15)
  mtext("Blue circles: A, all observations    Orange triangles: B, omit only E-46/R-3",
        side = 3, outer = TRUE, line = 0.5, cex = 0.95)
  mtext("Numeric entry order; intervals are for means, not individual observations or pairwise tests.",
        side = 1, outer = TRUE, line = 0.5, cex = 0.9)
  dev.off()
}
changes <- do.call(rbind, changes)
largest <- do.call(rbind, largest)
write.csv(changes, "outputs/entry_estimates_scenario_changes.csv", row.names = FALSE)
write.csv(largest, "outputs/entry_estimates_largest_sensitivity_changes.csv", row.names = FALSE)
cat("Validation checks:\n")
print(validation, row.names = FALSE)
cat("\nLargest absolute sensitivity changes (not a performance ranking):\n")
print(largest, row.names = FALSE, digits = 8)
cat("\nE-46 explicitly:\n")
print(all_estimates[all_estimates$entry == "E-46", 1:8], row.names = FALSE, digits = 8)
notes <- c(
  "# Entry estimates with equal block weighting", "",
  "Run `Rscript R/09_entry_estimates.R` (base R). Model definitions, factor levels, and data checks follow inspected scripts 07 and 08: plant_weight ~ replication + entry and sqrt(plant_weight) ~ replication + entry. A includes all 159 observations; B omits only E-46/R-3 as sensitivity, without declaring the zero erroneous. All four models are full rank.", "",
  "## Estimand and calculation", "",
  "For each entry, predict at R-1, R-2, and R-3 and average with weights 1/3 each. Construct those three design-matrix rows using the fitted model's terms, factor levels, and contrasts, then average them to obtain L. Estimate = L beta; SE = sqrt(L vcov(model) L'). The full covariance matrix, including off-diagonal covariances, is used. Averaging three individual prediction SEs would be incorrect.",
  "Equal weighting targets the same three-block average for every entry, even E-46 in B. Its R-3 response is model-estimated, not an imputed observation inserted into the dataset. An adjusted mean is a model-based average response, not an individual plant observation. Units and sampling details remain unconfirmed.",
  "Pointwise 95% confidence limits are estimate +/- qt(0.975, residual df) * SE. These are mean-response confidence intervals, not prediction intervals. Pointwise means each interval has its stated coverage under the assumptions when considered separately; the collection does not have simultaneous 95% coverage. Overlapping or nonoverlapping intervals are not a formal pairwise test.",
  "Square-root estimates and intervals stay on the modeled square-root scale. They are not squared or labelled arithmetic mean grain weight. recorded_yield and its unresolved conversion are not analyzed.", "",
  "## Verification", "",
  "All four models have 53 estimates in E-1 to E-53 numeric order. Estimates are checked against averaged predict() values. SEs are checked independently using the design matrix QR decomposition. Scenario A original-scale adjusted means agree with raw entry means within 1e-10, as expected for the balanced additive model.",
  capture.output(print(validation, row.names = FALSE)), "",
  "## Sensitivity changes", "",
  "The following table identifies the largest absolute A-to-B changes, not a breeding or performance ranking. Differences concern means on each respective response scale and are not tested for significance.",
  "E-46 is the only entry whose estimated mean changes beyond numerical tolerance (1e-10). Its original-scale increase is about 0.1252 and its square-root-scale increase about 0.6201. All other changes are floating-point noise (less than 2e-13), not meaningful differences between scenarios; their apparent order in the five-row table has no substantive interpretation. Balanced entries retain their three-block raw means under the additive model. Interval widths can still change when model variance and degrees of freedom change.",
  capture.output(print(largest, row.names = FALSE, digits = 8)), "",
  "E-46 estimates and intervals:",
  capture.output(print(all_estimates[all_estimates$entry == "E-46", 1:8], row.names = FALSE, digits = 8)), "",
  "## Interpretation and limits", "",
  "Intervals rely on the additive model, independent errors, common residual variance, and normal-error assumptions for exact t inference. Residual variance concerns from scripts 07 and 08 persist, so these are provisional model-based intervals, not assumption-free uncertainty estimates. Diagnostics cannot verify randomization or independence.",
  "The original-scale E-46 intervals extend below zero because the unconstrained linear-model t interval is symmetric. The limits are not clipped or interpreted as physically possible negative grain weights. Both scenario-comparison plots were inspected; A and B estimates largely coincide except for E-46, whose interval is wider in B. The square-root plot remains entirely on the modeled scale.",
  "These estimates concern the reported 2021 trial and these blocks, not broad adaptation across environments or years. They do not establish genetic superiority. Neither scenario nor response scale is selected based on significance. No pairwise tests, significance-group letters, or breeding shortlist is produced.", "",
  "## Outputs", "",
  "outputs/entry_estimates_all_models.csv contains all 212 estimates; four scenario/scale CSVs each contain 53 rows. entry_estimates_averaged_design_rows.csv records the L rows; entry_estimates_validation.csv records numerical checks. entry_estimates_scenario_changes.csv retains all changes in numeric order within each scale, and entry_estimates_largest_sensitivity_changes.csv reports the five largest absolute changes per scale. Two scenario_comparison PNGs show A and B estimates and pointwise intervals in numeric entry order.",
  "All existing scripts, model files, and source observations are preserved.", "", R.version.string
)
writeLines(notes, "notes/09_entry_estimates.md")
