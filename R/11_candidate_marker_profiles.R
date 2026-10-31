# Run from project folder: Rscript R/11_candidate_marker_profiles.R
# Base R descriptive mismatch comparison; no imputation or association tests.
read_output <- function(name) read.csv(file.path("outputs", name), stringsAsFactors = FALSE)
candidate_data <- read_output("candidate_sensitivity_comparison.csv")
mapping <- read_output("sample_id_mapping.csv")
samples <- read_output("kasp_sample_audit.csv")
markers <- read_output("kasp_marker_audit.csv")
calls <- read_output("kasp_original_and_analysis_calls.csv")
membership_columns <- c("A_original_candidate", "A_sqrt_candidate",
                        "B_original_candidate", "B_sqrt_candidate")
stopifnot(all(membership_columns %in% names(candidate_data)))
ten <- candidate_data[candidate_data$requested_size == 10, ]
stopifnot(nrow(ten) == 53, !anyDuplicated(ten$entry))
candidate_sets <- lapply(ten[membership_columns], function(x) ten$entry[x])
stopifnot(all(vapply(candidate_sets, function(x) setequal(x, candidate_sets[[1]]), TRUE)))
entries <- paste0("E-", 1:53)
candidates <- entries[entries %in% candidate_sets[[1]]]
stopifnot(length(candidates) == 10, "E-37" %in% candidates,
          !anyNA(mapping), !anyDuplicated(mapping$phenotype_entry),
          !anyDuplicated(mapping$kasp_id), !anyDuplicated(samples$kasp_id),
          setequal(mapping$phenotype_entry, entries), nrow(samples) == 53,
          setequal(mapping$kasp_id, samples$kasp_id), nrow(markers) == 502,
          !anyDuplicated(markers$excel_column),
          !anyDuplicated(calls[c("kasp_id", "excel_column")]))
# Explicit key matching, including sample-name and confirmed suffix checks.
mi <- match(candidates, mapping$phenotype_entry)
qi <- match(mapping$kasp_id[mi], samples$kasp_id)
stopifnot(!anyNA(mi), !anyNA(qi),
          all(mapping$kasp_id[mi] == paste0("NYUT_21_", sub("E-", "", candidates))),
          all(mapping$sample_name[mi] == samples$sample_name[qi]))
candidate_qc <- cbind(data.frame(phenotype_entry = candidates), samples[qi, ])
# These flags reproduce exploratory sample limits, not exclusions or validated QC.
for (limit in c(10, 20, 30)) {
  candidate_qc[[paste0("exceeds_exploratory_sample_missingness_", limit, "pct")]] <-
    candidate_qc$missing_calls > (limit / 100) * candidate_qc$total_markers
}
dir.create("outputs", showWarnings = FALSE)
dir.create("notes", showWarnings = FALSE)
write.csv(candidate_qc, "outputs/candidate_marker_profiles_sample_qc.csv", row.names = FALSE)

# Rebuild all 53 x 502 cells by keys. Revalidate canonical calls against originals.
i <- match(calls$kasp_id, samples$kasp_id)
j <- match(calls$excel_column, markers$excel_column)
stopifnot(!anyNA(i), !anyNA(j), nrow(calls) == 53 * 502)
original <- calls$original_call
missing <- is.na(original) | trimws(original) == "" | original %in% c("?", "Bad", "Uncallable")
valid <- grepl("^[ACGT]:[ACGT]$", original)
stopifnot(all(missing | valid))
canonical <- rep(NA_character_, length(original))
canonical[valid] <- vapply(strsplit(original[valid], ":", fixed = TRUE),
                          function(x) paste(sort(x), collapse = ":"), "")
stopifnot(identical(is.na(canonical), is.na(calls$analysis_call)),
          all(canonical[valid] == calls$analysis_call[valid]),
          all(calls$call_status == ifelse(missing, "missing", "valid")))
G <- matrix(NA_character_, 53, 502)
G[cbind(i, j)] <- canonical
missing_counts <- colSums(is.na(G))
alleles <- vapply(seq_len(502), function(k) {
  x <- G[, k]; length(unique(unlist(strsplit(x[!is.na(x)], ":", fixed = TRUE))))
}, integer(1))
stopifnot(all(missing_counts == markers$missing_calls),
          all(alleles == markers$number_observed_alleles),
          all(rowSums(is.na(G)) == samples$missing_calls))
candidate_rows <- match(mapping$kasp_id[mi], samples$kasp_id)
panels <- list()
pair_tables <- list()
panel_membership <- list()
pair_indices <- combn(seq_along(candidates), 2)
for (limit in c(10, 20, 30)) {
  panel <- paste0("M", limit)
  keep <- missing_counts <= (limit / 100) * 53 & alleles == 2
  panels[[panel]] <- which(keep)
  panel_membership[[panel]] <- data.frame(
    panel = panel, maximum_missingness = limit / 100,
    excel_column = markers$excel_column, marker_name = markers$marker_name,
    position = markers$position, original_sample_denominator = 53,
    missing_calls = missing_counts, original_missingness = missing_counts / 53,
    original_observed_alleles = alleles, in_panel = keep)
  results <- lapply(seq_len(ncol(pair_indices)), function(p) {
    a <- pair_indices[1, p]; b <- pair_indices[2, p]
    ga <- G[candidate_rows[a], keep]; gb <- G[candidate_rows[b], keep]
    jointly <- !is.na(ga) & !is.na(gb)
    n <- sum(jointly); mismatches <- sum(ga[jointly] != gb[jointly])
    data.frame(panel = panel, maximum_missingness = limit / 100,
      entry_1 = candidates[a], entry_2 = candidates[b],
      panel_markers = sum(keep), jointly_called_markers = n,
      fraction_panel_available = if (sum(keep) > 0) n / sum(keep) else NA_real_,
      mismatched_markers = mismatches,
      mismatch_proportion = if (n > 0) mismatches / n else NA_real_,
      no_joint_calls = n == 0)
  })
  pair_tables[[panel]] <- do.call(rbind, results)
}
pairs <- do.call(rbind, pair_tables)
# Relative review flags avoid asserting an unsupported reliability cutoff.
# Flag denominators or available fractions in the panel's lowest quartile.
for (panel in names(panels)) {
  rows <- which(pairs$panel == panel)
  z <- pairs[rows, ]
  qn <- unname(quantile(z$jointly_called_markers, 0.25))
  qf <- if (all(is.na(z$fraction_panel_available))) NA_real_ else
    unname(quantile(z$fraction_panel_available, 0.25, na.rm = TRUE))
  pairs$count_lower_quartile_reference[rows] <- qn
  pairs$fraction_lower_quartile_reference[rows] <- qf
  pairs$limited_shared_data_relative_flag[rows] <- z$no_joint_calls |
    z$jointly_called_markers <= qn |
    (!is.na(z$fraction_panel_available) & !is.na(qf) & z$fraction_panel_available <= qf)
}
write.csv(do.call(rbind, panel_membership), "outputs/candidate_marker_panel_membership.csv", row.names = FALSE)
write.csv(pairs, "outputs/candidate_marker_pairwise_profiles.csv", row.names = FALSE)

# Wide panel comparison, preserving the same 45 numeric-order pairs.
comparison <- pair_tables[[1]][c("entry_1", "entry_2")]
for (panel in names(panels)) {
  z <- pairs[pairs$panel == panel, ]
  stopifnot(identical(z$entry_1, comparison$entry_1), identical(z$entry_2, comparison$entry_2))
  for (field in c("mismatch_proportion", "jointly_called_markers", "fraction_panel_available", "limited_shared_data_relative_flag")) {
    comparison[[paste(panel, field, sep = "_")]] <- z[[field]]
  }
}
props <- comparison[paste0(names(panels), "_mismatch_proportion")]
comparison$panel_mismatch_range <- apply(props, 1, function(x) {
  if (sum(!is.na(x)) < 2) NA_real_ else diff(range(x, na.rm = TRUE))
})
write.csv(comparison, "outputs/candidate_marker_panel_comparison.csv", row.names = FALSE)

palette <- colorRampPalette(c("#F7FBFF", "#80B8D4", "#24638B", "#132D50"))(101)
for (panel in names(panels)) {
  z <- pairs[pairs$panel == panel, ]
  mismatch <- counts <- matrix(NA_real_, 10, 10, dimnames = list(candidates, candidates))
  for (p in seq_len(nrow(z))) {
    a <- match(z$entry_1[p], candidates); b <- match(z$entry_2[p], candidates)
    mismatch[a, b] <- mismatch[b, a] <- z$mismatch_proportion[p]
    counts[a, b] <- counts[b, a] <- z$jointly_called_markers[p]
  }
  # Diagonals are intentionally not pairwise comparisons.
  write.csv(data.frame(entry = candidates, counts, check.names = FALSE),
            paste0("outputs/candidate_marker_", panel, "_joint_counts.csv"), row.names = FALSE)
  write.csv(data.frame(entry = candidates, mismatch, check.names = FALSE),
            paste0("outputs/candidate_marker_", panel, "_mismatch_matrix.csv"), row.names = FALSE)
  png(paste0("outputs/candidate_marker_", panel, "_heatmap.png"), width = 1700, height = 1500, res = 160)
  par(mar = c(5, 5, 5, 7), xpd = FALSE)
  plot(NA, xlim = c(0.5, 10.5), ylim = c(0.5, 10.5), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", asp = 1,
       main = paste0(panel, ": observed canonical-genotype mismatch"))
  axis(1, at = 1:10, labels = candidates, las = 2)
  axis(2, at = 10:1, labels = candidates, las = 1)
  for (a in 1:10) for (b in 1:10) {
    value <- mismatch[a, b]
    color <- if (is.na(value)) "gray93" else palette[1 + round(100 * value)]
    rect(b - 0.5, 11 - a - 0.5, b + 0.5, 11 - a + 0.5, col = color, border = "white")
    text(b, 11 - a, if (a == b) "—" else if (is.na(value)) "NA" else sprintf("%.2f", value),
         cex = 0.8, col = if (!is.na(value) && value > 0.5) "white" else "#173449")
  }
  par(xpd = NA)
  for (k in 1:100) rect(10.9, 0.5 + (k - 1) / 10, 11.2, 0.5 + k / 10, col = palette[k], border = NA)
  text(11.35, c(0.5, 5.5, 10.5), c("0", "0.5", "1"), pos = 4, cex = 0.85)
  mtext(paste(length(panels[[panel]]), "biallelic panel markers; denominators vary by pair. Diagonal not compared."),
        side = 3, line = 0.5, cex = 0.85)
  dev.off()
}
summary <- do.call(rbind, lapply(names(panels), function(panel) {
  z <- pairs[pairs$panel == panel, ]
  data.frame(panel = panel, panel_markers = length(panels[[panel]]), pairs = nrow(z),
    min_joint_calls = min(z$jointly_called_markers), max_joint_calls = max(z$jointly_called_markers),
    min_fraction_available = min(z$fraction_panel_available),
    max_fraction_available = max(z$fraction_panel_available),
    min_mismatch = min(z$mismatch_proportion, na.rm = TRUE),
    max_mismatch = max(z$mismatch_proportion, na.rm = TRUE),
    relative_review_flags = sum(z$limited_shared_data_relative_flag),
    no_joint_calls = sum(z$no_joint_calls))
}))
write.csv(summary, "outputs/candidate_marker_panel_summary.csv", row.names = FALSE)
print(summary, row.names = FALSE)
cat("\nLargest across-panel mismatch changes (descriptive sensitivity):\n")
largest <- comparison[head(order(comparison$panel_mismatch_range, decreasing = TRUE), 5), ]
print(largest[c("entry_1", "entry_2", paste0(names(panels), "_mismatch_proportion"), "panel_mismatch_range")], row.names = FALSE)
notes <- c(
  "# Candidate marker profiles", "",
  "Run `Rscript R/11_candidate_marker_profiles.R` (base R). Existing CSV structures were inspected first. All four top-ten membership sets from script 10 agree; all ten candidates, including E-37, remain included. Mapping is validated by explicit IDs, confirmed suffixes, and sample names. Correspondence remains user-confirmed, not independently field-verified.", "",
  "## Definitions", "",
  "Panels M10, M20, and M30 retain markers with missing calls / 53 original samples <= 0.10, 0.20, or 0.30, respectively, and exactly two observed alleles across the same original 53 samples. These are exploratory panels. Missingness and allele counts are recalculated from the audited calls and checked against audit 04. Sample filtering is never applied. The membership CSV includes every original marker for every panel.",
  "Original genotype calls are read alongside their documented canonical representation and revalidated. Valid A/C/G/T pairs are alphabetically sorted: G:A and A:G are equivalent. A:A and A:G remain distinct genotypes. ?, Uncallable, Bad, and blanks remain missing. Unexpected strings stop the script. No calls are imputed.",
  "For each of 45 unordered candidate pairs and each panel, the denominator is the number of panel markers with valid calls in BOTH samples. Mismatch proportion = differing canonical genotypes / jointly called markers. No joint calls yields NA. The available fraction is jointly called markers / all panel markers, not / 502. Different pairs and panels can use different marker subsets.",
  "Limited-data review flags identify denominators or available fractions at or below the panel-specific lower quartile among the 45 pairs, plus any pair with no joint calls. The quartile reference values and continuous counts/fractions are exported. This is an explicitly relative descriptive flag, NOT a validated reliability cutoff; unflagged comparisons are not certified reliable.",
  "Candidate QC metrics and flags for exceeding the earlier exploratory 10%, 20%, and 30% sample-missingness limits are preserved in a separate CSV. No candidate is removed because of QC. Missingness uses 502 original markers; heterozygosity uses valid calls as denominator.", "",
  "## Actual results", "", capture.output(print(summary, row.names = FALSE)), "",
  "Largest across-panel changes (not a ranking of crosses):", "",
  capture.output(print(largest[c("entry_1", "entry_2", paste0(names(panels), "_mismatch_proportion"), "panel_mismatch_range")], row.names = FALSE)), "",
  "## Interpretation and outputs", "",
  "All three heatmaps were visually inspected. Broad patterns persist across panels, with E-11 showing higher mismatch against other candidates than most other displayed pairs. E-42/E-45 has zero observed mismatch across 363, 430, and 443 jointly called markers in M10, M20, and M30. This does not establish genome-wide identity, duplicate sample status, or pedigree. The largest across-panel range is about 0.0410 for E-3/E-45 (4.10 percentage points).",
  "Each panel has 12 relative lower-quartile review flags, including all nine pairs involving E-37. E-37 remains retained, with its original QC metrics and exploratory flags exported. A lower available fraction means a larger part of that panel could not contribute to the pair comparison; a zero mismatch is agreement among observed calls only.",
  "These quantities describe observed genotype mismatch in a targeted marker panel. They are not genome-wide genetic distances, causal marker effects, favorable alleles, or evidence that a cross will produce heterosis. Panel selection and missing data affect which loci contribute, so variation across panels is descriptive sensitivity evidence. No crosses are recommended from mismatch alone.",
  "All three heatmaps use the same fixed numeric entry order and the same 0–1 color scale, without clustering. Cell labels show mismatch proportions rounded to two decimals. Diagonals are not pair comparisons. Joint-call matrices and long reports provide exact denominators and available fractions.",
  "Outputs under outputs/: candidate_marker_panel_membership.csv; candidate_marker_pairwise_profiles.csv (135 rows); candidate_marker_panel_comparison.csv (45 rows); candidate_marker_panel_summary.csv; candidate_marker_profiles_sample_qc.csv; and M10/M20/M30 mismatch matrices, joint-count matrices, and heatmap PNGs with the candidate_marker_ prefix.",
  "Source observations and prior outputs are preserved. No candidates are removed, calls imputed, associations tested, or crosses recommended.", "", R.version.string
)
writeLines(notes, "notes/11_candidate_marker_profiles.md")
