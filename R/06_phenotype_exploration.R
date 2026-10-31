# Run from the project folder: Rscript R/06_phenotype_exploration.R
# Base R descriptive exploration only: observations remain unchanged.
phenotype <- read.csv("outputs/phenotype_2021_imported.csv",
                      stringsAsFactors = FALSE, check.names = FALSE)
measurements <- c("plant_height", "tillers", "panicle_length", "filled_grains",
                  "unfilled_grains", "grain_weight_100", "plant_weight", "recorded_yield")
expected_columns <- c("entry", measurements, "entry_original", "replication")
if (!identical(names(phenotype), expected_columns)) stop("Unexpected phenotype columns.")
if (!all(vapply(phenotype[measurements], is.numeric, logical(1)))) {
  stop("Every measurement column must be numeric; no automatic conversion applied.")
}
if (any(vapply(phenotype[measurements], function(x) any(!is.finite(x) & !is.na(x)), logical(1)))) {
  stop("An infinite measurement requires inspection before exploration.")
}
replications <- c("R-1", "R-2", "R-3")
stopifnot(!anyNA(phenotype[c("entry", "replication")]),
          all(nzchar(phenotype$entry)),
          setequal(unique(phenotype$replication), replications),
          !anyDuplicated(phenotype[c("entry", "replication")]))
highlight <- phenotype$entry == "E-46" & phenotype$replication == "R-3"
stopifnot(sum(highlight) == 1, phenotype$plant_weight[highlight] == 0)
cat("Validated columns and numeric measurements:", nrow(phenotype), "records.\n")

# Grouping means selecting rows belonging to the same replication.
# For each group and measurement, NA values are counted then omitted from
# descriptive calculations. Recorded zeros are included in all statistics.
summaries <- list()
for (replication_name in replications) {
  group <- phenotype[phenotype$replication == replication_name, ]
  for (measurement in measurements) {
    x <- group[[measurement]]
    observed <- x[!is.na(x)]
    n <- length(observed)
    summaries[[paste(replication_name, measurement)]] <- data.frame(
      replication = replication_name, measurement = measurement,
      number_observed = n, missing_count = sum(is.na(x)),
      zero_count = sum(x == 0, na.rm = TRUE),
      mean = if (n > 0) mean(observed) else NA_real_,
      standard_deviation = if (n > 1) sd(observed) else NA_real_,
      minimum = if (n > 0) min(observed) else NA_real_,
      median = if (n > 0) median(observed) else NA_real_,
      maximum = if (n > 0) max(observed) else NA_real_
    )
  }
}
summary_by_replication <- do.call(rbind, summaries)
rownames(summary_by_replication) <- NULL
grain_summary <- summary_by_replication[
  summary_by_replication$measurement == "plant_weight", ]
dir.create("outputs", showWarnings = FALSE)
dir.create("notes", showWarnings = FALSE)
write.csv(summary_by_replication, "outputs/phenotype_summary_by_replication.csv", row.names = FALSE)
write.csv(grain_summary, "outputs/grain_weight_per_plant_summary_by_replication.csv", row.names = FALSE)
cat("\nGrain weight per plant (unit unconfirmed):\n")
print(grain_summary, row.names = FALSE, digits = 6)

# The user confirms plant_weight means grain weight per plant; units are
# unconfirmed. recorded_yield remains a separate measurement column.
x_replication <- match(phenotype$replication, replications)
observed_weights <- phenotype$plant_weight[!is.na(phenotype$plant_weight)]
plot_range <- range(c(0, observed_weights))
span <- max(diff(plot_range), 1)
plot_range <- plot_range + c(-0.04, 0.12) * span
accent <- "#BD3D25"
png("outputs/grain_weight_per_plant_by_replication.png", width = 1600, height = 1100, res = 150)
par(mar = c(5, 5, 5, 2))
boxplot(lapply(replications, function(r) phenotype$plant_weight[phenotype$replication == r]),
        names = replications, outline = FALSE, col = "#E5EDF0", border = "#667B84",
        ylim = plot_range, xlab = "Replication",
        ylab = "Grain weight per plant (unit unconfirmed)",
        main = "Grain weight per plant by replication")
# Jitter spreads points horizontally for visibility, never changing their y values.
set.seed(2021)
points(jitter(x_replication, amount = 0.15), phenotype$plant_weight,
       pch = 16, col = adjustcolor("#286478", alpha.f = 0.6), cex = 0.8)
points(x_replication[highlight], phenotype$plant_weight[highlight],
       pch = 21, bg = accent, col = "white", cex = 1.6)
text(3, 0, "E-46 / R-3: recorded zero", pos = 2, offset = 0.8, col = accent, cex = 0.85)
mtext("Boxes summarize observations; points show every record. Zero retained unchanged.",
      side = 3, line = 0.5, cex = 0.85)
dev.off()

# Profiles are assembled by entry AND replication keys, not input row order.
entries <- unique(phenotype$entry)
profiles <- matrix(NA_real_, nrow = length(entries), ncol = 3,
                   dimnames = list(entries, replications))
profiles[cbind(match(phenotype$entry, entries), x_replication)] <- phenotype$plant_weight
png("outputs/grain_weight_per_plant_profiles.png", width = 1600, height = 1100, res = 150)
par(mar = c(5, 5, 5, 2))
plot(1:3, rep(NA_real_, 3), type = "n", xaxt = "n", ylim = plot_range,
     xlim = c(0.9, 3.2), xlab = "Replication",
     ylab = "Grain weight per plant (unit unconfirmed)",
     main = "Grain weight per plant: entry profiles")
axis(1, at = 1:3, labels = replications)
for (entry in setdiff(entries, "E-46")) {
  lines(1:3, profiles[entry, ], col = adjustcolor("#667B84", alpha.f = 0.27), lwd = 1)
}
lines(1:3, profiles["E-46", ], type = "b", pch = 16, col = accent, lwd = 2.5)
text(3, profiles["E-46", 3], "E-46: recorded zero", pos = 2, offset = 0.8,
     col = accent, cex = 0.85)
legend("topright", legend = c("E-46", "Other entries"),
       col = c(accent, "#9AA7AC"), lwd = c(2.5, 1), bty = "n")
mtext("Lines connect the same entry across replications; they do not represent time trends.",
      side = 3, line = 0.5, cex = 0.85)
dev.off()

table_rows <- sprintf("| %s | %d | %d | %d | %.3f | %.3f | %.3f | %.3f | %.3f |",
  grain_summary$replication, grain_summary$number_observed, grain_summary$missing_count,
  grain_summary$zero_count, grain_summary$mean, grain_summary$standard_deviation,
  grain_summary$minimum, grain_summary$median, grain_summary$maximum)
e46_text <- paste(sprintf("%s = %.3f", replications, profiles["E-46", ]), collapse = "; ")
notes <- c(
  "# Phenotype descriptive exploration", "",
  "Run `Rscript R/06_phenotype_exploration.R` from the project folder. Only base R is used. The expected 11 columns and eight numeric measurement columns are validated before summarizing.", "",
  "## Meaning and unresolved questions", "",
  "The user confirms that plant_weight means grain weight per plant. Its physical unit remains unconfirmed, so both plots label it 'unit unconfirmed'. Original column names and data are preserved.",
  "E-46/R-3 has a recorded zero. It is retained in summaries and plots without interpreting it as confirmed crop failure or as missing data. Its biological meaning requires clarification.",
  "recorded_yield is summarized separately, never substituted for grain weight per plant. The factor 135.85 is an unverified numerical relationship, not a validated conversion to yield per hectare. Audit 02 found factor discrepancies for E-31/R-1 (+1000), E-39/R-3 (+1000), and E-40/R-3 (-3000), rounded. Combined differs from replication recorded_yield for E-42/R-2 (-2000) and E-46/R-3 (+130), expressed as Combined minus replication. See notes/02_yield_consistency.md. No correction is made here.", "",
  "## Grain weight per plant summaries", "",
  "| Replication | Observed | Missing | Zeros | Mean | SD | Minimum | Median | Maximum |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|", table_rows, "",
  "## How to read the plots", "",
  "The replication plot shows every individual observation, with horizontal jitter to reduce overlap. The vertical values are unchanged. Boxes show the median and middle half of the data; whiskers extend to the most extreme observations within 1.5 times the interquartile range. All observations, including those beyond whiskers, appear as points. The E-46/R-3 zero is highlighted.",
  paste("The profile plot highlights E-46:", e46_text),
  "Other entries are subdued but all are shown. Connected profiles reveal within-entry variation across the three replications. Replication is a grouping label, not a time axis. Crossing lines show that entries do not follow identical profiles.",
  "Differences between the displayed means, spreads, and profiles are descriptive observations only. These plots establish neither statistically significant replication differences nor genetic superiority. No lines are ranked.", "",
  "## Functions and denominators", "",
  "Grouping uses phenotype[phenotype$replication == replication_name, ] to select one replication's rows. Each measurement is then summarized separately within that group.",
  "mean() adds the observed values and divides by their count. sd() measures spread around the mean, using the sample standard deviation with denominator n - 1. median() returns the middle sorted value (or the average of the two middle values for an even count).",
  "Missing values are counted separately and omitted from these calculations; zeros are observed values and remain included. SD is NA with fewer than two observations; other summaries are NA when there are no observations. The output records number_observed, missing_count, and zero_count explicitly.", "",
  "## Outputs", "",
  "- outputs/phenotype_summary_by_replication.csv: all eight measurements, each summarized in all three replications.",
  "- outputs/grain_weight_per_plant_summary_by_replication.csv: the plant_weight subset.",
  "- outputs/grain_weight_per_plant_by_replication.png: distributions and individual observations.",
  "- outputs/grain_weight_per_plant_profiles.png: one connected profile per entry.",
  "No observations were removed or replaced, outliers filtered, or statistical models fitted. All previous files remain unchanged."
)
writeLines(notes, "notes/06_phenotype_exploration.md")
