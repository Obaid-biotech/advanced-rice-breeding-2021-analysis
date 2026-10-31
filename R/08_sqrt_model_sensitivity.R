# Run from the project folder: Rscript R/08_sqrt_model_sensitivity.R
# Base R only. Previous scripts, outputs, and source observations are preserved.
data <- read.csv("outputs/phenotype_2021_imported.csv", stringsAsFactors = FALSE)
stopifnot(all(c("entry", "replication", "plant_weight") %in% names(data)),
          is.numeric(data$plant_weight), nrow(data) == 159,
          !anyNA(data[c("entry", "replication", "plant_weight")]),
          all(is.finite(data$plant_weight)),
          !anyDuplicated(data[c("entry", "replication")]))
# Check before sqrt(): zero is allowed, negative values are not.
if (any(data$plant_weight < 0)) stop("Negative responses: square root is undefined here.")
stopifnot(setequal(data$entry, paste0("E-", 1:53)),
          setequal(data$replication, c("R-1", "R-2", "R-3")),
          all(table(data$entry, data$replication) == 1))
data$entry <- factor(data$entry, levels = paste0("E-", 1:53))
data$replication <- factor(data$replication, levels = c("R-1", "R-2", "R-3"))
omit <- data$entry == "E-46" & data$replication == "R-3"
stopifnot(sum(omit) == 1, data$plant_weight[omit] == 0)
scenarios <- list(A_all_recorded = data, B_omit_E46_R3 = data[!omit, ])
dir.create("outputs", showWarnings = FALSE)
dir.create("notes", showWarnings = FALSE)
summaries <- list()
tests <- list()
diagnostic_tables <- list()
for (scenario in names(scenarios)) {
  d <- scenarios[[scenario]]
  models <- list(
    original = lm(plant_weight ~ replication + entry, data = d, na.action = na.fail),
    sqrt = lm(sqrt(plant_weight) ~ replication + entry, data = d, na.action = na.fail)
  )
  # Save only new model copies; nothing from script 07 is overwritten.
  saveRDS(models, paste0("outputs/sqrt_sensitivity_", scenario, "_models.rds"))
  png(paste0("outputs/sqrt_sensitivity_", scenario, "_diagnostics.png"),
      width = 2600, height = 2250, res = 180)
  # Rows are diagnostic types; columns are response scales.
  par(mfrow = c(3, 2), mar = c(4.5, 4.5, 3, 1.5), oma = c(0, 0, 3, 0))
  diagnostic_data <- list()
  for (scale in names(models)) {
    model <- models[[scale]]
    key <- paste(scenario, scale, sep = "_")
    rank_ok <- model$rank == ncol(model.matrix(model))
    stopifnot(rank_ok, df.residual(model) > 0)
    summaries[[key]] <- data.frame(
      scenario = scenario, response_scale = scale, observations = nobs(model),
      model_rank = model$rank, design_columns = ncol(model.matrix(model)),
      full_rank = rank_ok, residual_df = df.residual(model),
      residual_sd_on_response_scale = sigma(model),
      negative_fitted_values = sum(fitted(model) < 0)
    )
    partial <- as.data.frame(drop1(model, test = "F"))
    # Omit ordinary AIC from cross-scale reports: responses differ.
    partial$AIC <- NULL
    partial$term <- rownames(partial)
    rownames(partial) <- NULL
    partial$scenario <- scenario
    partial$response_scale <- scale
    tests[[key]] <- partial
    diagnostic <- data.frame(
      scenario = scenario, response_scale = scale, entry = d$entry,
      replication = d$replication, observed_original = d$plant_weight,
      observed_on_response_scale = if (scale == "sqrt") sqrt(d$plant_weight) else d$plant_weight,
      fitted = fitted(model), residual = residuals(model),
      standardized_residual = rstandard(model)
    )
    # Label the three largest absolute standardized residuals in this fit,
    # plus E-46/R-3 if present. This changes labels only, never fitted data.
    label_index <- union(head(order(abs(diagnostic$standardized_residual), decreasing = TRUE), 3),
                         which(d$entry == "E-46" & d$replication == "R-3"))
    diagnostic$labelled <- seq_len(nrow(d)) %in% label_index
    diagnostic_data[[scale]] <- diagnostic
    diagnostic_tables[[key]] <- diagnostic
  }
  for (kind in c("residual_vs_fitted", "normal_qq", "scale_location")) {
    for (scale in names(models)) {
      z <- diagnostic_data[[scale]]
      labs <- paste(z$entry, z$replication, sep = "/")
      cols <- ifelse(z$entry == "E-46" & z$replication == "R-3", "#BE3928", "#31687B")
      units <- if (scale == "sqrt") "sqrt(grain weight per plant)" else "Grain weight per plant"
      title <- paste(if (scale == "sqrt") "Square-root response" else "Original response",
                     gsub("_", " ", kind), sep = ": ")
      if (kind == "normal_qq") {
        ord <- order(z$standardized_residual)
        x <- qnorm(ppoints(nrow(z)))
        y <- z$standardized_residual[ord]
        yr <- range(y); yr <- yr + c(-0.08, 0.15) * diff(yr)
        plot(x, y, pch = 16, col = cols[ord], main = title, ylim = yr,
             xlab = "Theoretical normal quantiles", ylab = "Standardized residuals")
        qqline(z$standardized_residual, col = "#B47A22", lwd = 2)
        pick <- which(z$labelled[ord])
        text(x[pick], y[pick], labs[ord][pick], pos = 3, cex = 0.7,
             col = cols[ord][pick], xpd = NA)
      } else {
        x <- z$fitted
        y <- if (kind == "scale_location") sqrt(abs(z$standardized_residual)) else z$residual
        yr <- range(y); yr <- yr + c(-0.08, 0.15) * diff(yr)
        plot(x, y, pch = 16, col = cols, main = title, ylim = yr,
             xlab = paste("Fitted", units, "(unit unconfirmed)"),
             ylab = if (kind == "scale_location") "sqrt(|standardized residual|)" else "Residual on fitted response scale")
        if (kind == "residual_vs_fitted") abline(h = 0, lty = 2, col = "gray60")
        lines(lowess(x, y), col = "#B47A22", lwd = 2)
        pick <- which(z$labelled)
        text(x[pick], y[pick], labs[pick], pos = 3, cex = 0.7, col = cols[pick], xpd = NA)
      }
    }
  }
  mtext(paste(scenario, "— diagnostic comparison (not back-transformed residuals)"),
        outer = TRUE, side = 3, line = 1, cex = 1.2)
  dev.off()
}
summary_table <- do.call(rbind, summaries)
partial_tests <- do.call(rbind, tests)
rownames(summary_table) <- rownames(partial_tests) <- NULL
write.csv(summary_table, "outputs/sqrt_sensitivity_model_summary.csv", row.names = FALSE)
write.csv(partial_tests, "outputs/sqrt_sensitivity_partial_tests.csv", row.names = FALSE)
write.csv(do.call(rbind, diagnostic_tables), "outputs/sqrt_sensitivity_diagnostics.csv", row.names = FALSE)
print(summary_table, row.names = FALSE)
print(partial_tests[partial_tests$term != "<none>", ], row.names = FALSE)
notes <- c(
  "# Square-root response sensitivity", "",
  "Run `Rscript R/08_sqrt_model_sensitivity.R` (base R only). All 159 responses were checked to be finite and nonnegative before square roots. The recorded zero is permitted and unchanged.", "",
  "## Models and interpretation", "",
  "Original: plant_weight ~ replication + entry. Transformed: sqrt(plant_weight) ~ replication + entry. Both predictors are factors. These remain provisional additive fixed-effects analyses of the reported 2021 RCBD, with unconfirmed units and sampling details.",
  "A includes all 159 observations. B omits only E-46/R-3, explicitly as a sensitivity analysis, without declaring it erroneous. Original fits are refitted for comparison without replacing script 07 outputs.",
  "Transformation changes the response scale and interpretation: the transformed model describes mean square-root grain weight, not arithmetic mean grain weight. Squaring a fitted square-root mean does not directly recover the arithmetic mean grain weight. For a transformed variable Z, E(Z^2) = [E(Z)]^2 + Var(Z); simply squaring omits the variance component. No back-transformed mean estimates are produced here.",
  "Residual SD is in original units for the original fit and square-root units for the transformed fit; a smaller numerical SD across these scales is not evidence of a better model. Ordinary lm AIC values must not be compared directly across different response scales without accounting for the transformation of the response density. Cross-scale AIC is deliberately omitted from the exported comparison.",
  "Model choice must not depend on obtaining a smaller p-value. The scientific response, design, sampling, and diagnostics matter. No preferred model is selected.", "",
  "## Results", "",
  capture.output(print(summary_table, row.names = FALSE)), "",
  capture.output(print(partial_tests[partial_tests$term != "<none>", ], row.names = FALSE)), "",
  "All four models are full rank. Partial F-tests use drop1(model, test = 'F'), testing each term while retaining the other. This avoids sequential term-order effects in unbalanced scenario B. The entry null is equal entry effects on the respective model's response scale after accounting for blocks; the response scales therefore have different null interpretations.", "",
  "## Diagnostic layout", "",
  "Each scenario's PNG has original-scale model diagnostics on the left and square-root model diagnostics on the right. Rows show residual-versus-fitted, normal Q-Q, and scale-location plots. Residuals are on each model's own response scale, not back-transformed. Only the three largest absolute standardized residuals per fit plus E-46/R-3, when present, are labelled. Smooth curves are descriptive guides.",
  "Diagnostics cannot verify randomization or independence. Field records, layout, experimental-unit definitions, and sampling details are still needed. No assumptions are declared solved merely because a transformation was fitted.", "",
  "## Visual inspection findings", "",
  "Both scenario comparison PNGs were inspected. The square-root fits reduce the pronounced upward scale-location pattern seen in the original fits, especially across the central range, but some increase remains at higher fitted values. In scenario A the scale-location curve rises again at very low fitted values, where there are few observations; E-46/R-3 becomes a relatively large negative standardized residual rather than the near-zero residual seen on the original scale. Scenario B has a flatter low-to-middle fitted range after transformation, with remaining upward spread at the high end.",
  "The square-root Q-Q plots look closer to the reference overall than the original fits, although departures and extreme residuals remain. All four fits have nonnegative fitted values in this dataset. The transformation changes which observations look unusual and has not established that constant variance or normal errors hold. Neither transformation nor omitting the zero verifies randomization, independence, or the reported experimental design.",
  "The square-root entry F-statistics are approximately 5.885 (A) and 4.812 (B); replication F-statistics are 9.912 and 8.822. Their smaller p-values than the original fits are not a reason to select the transformation. The source zero is unchanged, and no further records are omitted.", "",
  "## Outputs and preservation", "",
  "outputs/sqrt_sensitivity_model_summary.csv, sqrt_sensitivity_partial_tests.csv, and sqrt_sensitivity_diagnostics.csv contain all four fits' summaries, tests, and observation diagnostics. Each scenario also has a diagnostic comparison PNG and a new models.rds containing its two fitted models.",
  "No additional observations were deleted, lines ranked, pairwise tests run, or source observations modified. Previous scripts and outputs are preserved.", "",
  R.version.string
)
writeLines(notes, "notes/08_sqrt_model_sensitivity.md")
