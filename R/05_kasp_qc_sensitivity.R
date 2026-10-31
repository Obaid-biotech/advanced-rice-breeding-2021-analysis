# Run from the project folder: Rscript R/05_kasp_qc_sensitivity.R
# Base R only. These scenarios explore sensitivity, not universal QC rules.
samples <- read.csv("outputs/kasp_sample_audit.csv", stringsAsFactors = FALSE)
markers <- read.csv("outputs/kasp_marker_audit.csv", stringsAsFactors = FALSE)
calls <- read.csv("outputs/kasp_original_and_analysis_calls.csv",
                  stringsAsFactors = FALSE)

# Reconstruct matrices by explicit sample and Excel-column keys, not row order.
stopifnot(!anyNA(samples$kasp_id), !anyDuplicated(samples$kasp_id),
          !anyNA(markers$excel_column), !anyDuplicated(markers$excel_column),
          !anyDuplicated(calls[c("kasp_id", "excel_column")]))
i <- match(calls$kasp_id, samples$kasp_id)
j <- match(calls$excel_column, markers$excel_column)
stopifnot(!anyNA(i), !anyNA(j), nrow(calls) == nrow(samples) * nrow(markers),
          all(calls$call_status %in% c("valid", "missing")))
# Unexpected calls need an explicit decision before scenario calculations.
# Audit 04 found none; this check prevents silently treating them as missing.
missing <- matrix(FALSE, nrow(samples), nrow(markers))
analysis_calls <- matrix(NA_character_, nrow(samples), nrow(markers))
missing[cbind(i, j)] <- calls$call_status == "missing"
analysis_calls[cbind(i, j)] <- calls$analysis_call
stopifnot(all(is.na(analysis_calls) == missing),
          all(grepl("^[ACGT]:[ACGT]$", analysis_calls[!missing])),
          all(samples$total_markers == ncol(missing)),
          all(rowSums(missing) == samples$missing_calls),
          all(abs(rowMeans(missing) - samples$missing_call_proportion) < 1e-12),
          all(colSums(missing) == markers$missing_calls))

limits <- c(0.10, 0.20, 0.30)
results <- list()
details <- list()
for (sample_limit in limits) {
  # Start afresh from original missingness over ALL 502 markers every time.
  # Compare integer counts with the limit times denominator to include equality.
  keep_sample <- samples$missing_calls <= sample_limit * ncol(missing)
  retained_n <- sum(keep_sample)
  excluded <- samples$kasp_id[!keep_sample]
  for (marker_limit in limits) {
    scenario <- paste0("S", round(100 * sample_limit), "_M", round(100 * marker_limit))
    marker_missing <- rep(NA_real_, ncol(missing))
    allele_n <- rep(NA_integer_, ncol(missing))
    keep_marker <- rep(FALSE, ncol(missing))
    if (retained_n > 0) {
      # Denominator is the number of eligible samples, not the original 53.
      missing_counts <- colSums(missing[keep_sample, , drop = FALSE])
      marker_missing <- missing_counts / retained_n
      keep_marker <- missing_counts <= marker_limit * retained_n
      # Recount distinct allele letters after sample exclusion, including
      # markers that fail the marker limit so every marker has an audit row.
      allele_n <- vapply(seq_len(ncol(missing)), function(k) {
        x <- analysis_calls[keep_sample, k]
        length(unique(unlist(strsplit(x[!is.na(x)], ":", fixed = TRUE))))
      }, integer(1))
    }
    results[[scenario]] <- data.frame(
      scenario = scenario, max_sample_missingness = sample_limit,
      max_marker_missingness = marker_limit,
      original_samples = nrow(samples), original_markers = ncol(missing),
      sample_missingness_denominator = ncol(missing),
      marker_missingness_denominator = retained_n,
      retained_samples = retained_n, retained_markers = sum(keep_marker),
      retained_polymorphic_biallelic_markers = sum(keep_marker & allele_n == 2, na.rm = TRUE),
      no_retained_samples = retained_n == 0,
      excluded_sample_ids = paste(excluded, collapse = ";")
    )
    details[[scenario]] <- data.frame(
      scenario = scenario, excel_column = markers$excel_column,
      marker_name = markers$marker_name, eligible_samples = retained_n,
      marker_missingness = marker_missing, observed_allele_count = allele_n,
      retained_by_missingness = keep_marker
    )
  }
}
comparison <- do.call(rbind, results)
rownames(comparison) <- NULL
write.csv(comparison, "outputs/kasp_qc_sensitivity.csv", row.names = FALSE)
write.csv(do.call(rbind, details), "outputs/kasp_qc_sensitivity_marker_details.csv",
          row.names = FALSE)
display <- comparison[c("scenario", "retained_samples", "retained_markers",
                        "retained_polymorphic_biallelic_markers", "no_retained_samples")]
print(display, row.names = FALSE)
cat("\nExcluded samples by original sample-missingness limit:\n")
for (limit in limits) {
  row <- comparison[comparison$max_sample_missingness == limit, ][1, ]
  cat(sprintf("%d%%: %s\n", round(limit * 100), row$excluded_sample_ids))
}

# Two panels use separate count axes: samples and markers have different scales.
png("outputs/kasp_qc_sensitivity.png", width = 1800, height = 1050, res = 150)
par(mfrow = c(1, 2), mar = c(9, 4, 4, 1), oma = c(3, 0, 3, 0))
colors <- rep(c("#287A8E", "#D38A26", "#7855A1"), each = 3)
draw_counts <- function(values, title, maximum) {
  x <- barplot(values, names.arg = comparison$scenario, las = 2,
               col = colors, border = NA, ylim = c(0, maximum * 1.13),
               main = title, ylab = "Count retained", cex.names = 0.9)
  text(x, values, labels = values, pos = 3, cex = 0.85)
  abline(h = maximum, lty = 3, col = "gray45")
}
draw_counts(comparison$retained_samples, "Samples (original total: 53)", nrow(samples))
draw_counts(comparison$retained_polymorphic_biallelic_markers,
            "Polymorphic biallelic markers", nrow(markers))
mtext("Exploratory KASP QC sensitivity", outer = TRUE, side = 3, cex = 1.3)
mtext("S = maximum original sample missingness (%); M = maximum marker missingness (%)",
      outer = TRUE, side = 1, cex = 0.85)
dev.off()

# Generate readable notes with the actual results from this run.
table_rows <- sprintf("| %d%% | %d%% | %d | %d | %d |",
  round(comparison$max_sample_missingness * 100),
  round(comparison$max_marker_missingness * 100), comparison$retained_samples,
  comparison$retained_markers, comparison$retained_polymorphic_biallelic_markers)
excluded_notes <- vapply(limits, function(limit) {
  row <- comparison[comparison$max_sample_missingness == limit, ][1, ]
  sprintf("- Sample limit %d%%: %s.", round(limit * 100), row$excluded_sample_ids)
}, "")
notes <- c(
  "# KASP QC sensitivity", "",
  "Run `Rscript R/05_kasp_qc_sensitivity.R`. This script uses base R and the original audit CSVs from step 04. All input files remain unchanged.", "",
  "These nine combinations are exploratory sensitivity scenarios, not validated universal thresholds. No final QC scenario is selected.", "",
  "## Filtering order and denominators", "",
  "1. Select eligible samples using their ORIGINAL missing-call count divided by all 502 markers. Maximum limits are inclusive (<=).",
  "2. For each of the 502 markers, recalculate missing-call proportion using only those eligible samples. Its denominator is the retained sample count shown below, not 53.",
  "3. Retain markers at or below the scenario's marker-missingness limit.",
  "4. Recalculate observed allele counts from valid calls among eligible samples. Count retained markers with exactly two observed alleles as polymorphic biallelic markers.", "",
  "Each scenario starts from the original audit. Sample eligibility is not recalculated after marker selection, and there is no iterative filtering or imputation. Missing calls retain the audit definition: ?, Uncallable, Bad, or blank. Unexpected strings cause the script to stop rather than silently changing that definition.",
  "If a scenario has no eligible samples, marker missingness and allele counts are undefined (NA in the detail report), retained counts are zero, and no_retained_samples is TRUE.", "",
  "## Actual comparison", "",
  "| Max sample missingness | Max marker missingness | Samples | Markers | Polymorphic biallelic markers |",
  "|---|---|---:|---:|---:|", table_rows, "",
  paste("Scenarios with no retained samples:", sum(comparison$no_retained_samples)), "",
  "Excluded sample IDs (identical for all three marker limits within each sample limit):", "",
  excluded_notes, "",
  "## Interpretation", "",
  "Tighter sample limits keep fewer samples. Marker missingness can change when samples are removed, so it must be recomputed before applying the marker limit. Within a fixed sample limit, relaxing the marker limit can retain more markers.",
  "Observed polymorphism also depends on which samples remain. If the only sample carrying an alternative allele is excluded, a formerly biallelic marker can become monomorphic. Original full-sample allele counts therefore cannot be reused. Counts here measure observed alleles, not proof that an allele is absent from the wider population.",
  "Retained marker counts include monomorphic markers that pass missingness; polymorphic biallelic counts are a separate descriptive subset. Neither count determines the best scenario on its own.", "",
  "## Files", "",
  "- outputs/kasp_qc_sensitivity.csv: all nine scenarios, denominators, counts, no-sample flags, and excluded IDs.",
  "- outputs/kasp_qc_sensitivity_marker_details.csv: recalculated missingness and allele counts for every marker in every scenario.",
  "- outputs/kasp_qc_sensitivity.png: retained samples and polymorphic biallelic markers in separate panels. Dotted lines show original totals (53 samples and 502 markers).",
  "No phenotype records were excluded, no marker calls were joined to phenotypes, and no marker–trait associations were tested.", "",
  R.version.string
)
writeLines(notes, "notes/05_kasp_qc_sensitivity.md")
