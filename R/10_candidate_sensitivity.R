# Run from the project folder: Rscript R/10_candidate_sensitivity.R
# Base R exploratory prioritization, not a test of superior lines.
est <- read.csv("outputs/entry_estimates_all_models.csv", stringsAsFactors = FALSE)
mapping <- read.csv("outputs/sample_id_mapping.csv", stringsAsFactors = FALSE)
qc <- read.csv("outputs/kasp_sample_audit.csv", stringsAsFactors = FALSE)
entries <- paste0("E-", 1:53)
scenarios <- c("A_all_recorded", "B_omit_E46_R3")
scales <- c("original", "sqrt")
tolerance <- 1e-10  # Absolute tolerance on each model's own response scale.
stopifnot(all(c("scenario", "response_scale", "entry", "estimate", "lower_95", "upper_95") %in% names(est)),
          nrow(est) == 212, !anyNA(est),
          !anyDuplicated(est[c("scenario", "response_scale", "entry")]),
          setequal(est$scenario, scenarios), setequal(est$response_scale, scales),
          all(vapply(est[c("estimate", "lower_95", "upper_95")], is.numeric, logical(1))),
          all(is.finite(as.matrix(est[c("estimate", "lower_95", "upper_95")]))))
for (s in scenarios) for (scale in scales) {
  stopifnot(setequal(est$entry[est$scenario == s & est$response_scale == scale], entries))
}
# Validate one-to-one mapping and QC joins, never assuming matching row orders.
stopifnot(all(c("phenotype_entry", "kasp_id", "sample_name") %in% names(mapping)),
          all(c("kasp_id", "sample_name", "missing_call_proportion", "heterozygous_call_proportion") %in% names(qc)),
          !anyNA(mapping), !anyNA(qc$kasp_id),
          !anyDuplicated(mapping$phenotype_entry), !anyDuplicated(mapping$kasp_id),
          !anyDuplicated(qc$kasp_id), setequal(mapping$phenotype_entry, entries),
          setequal(mapping$kasp_id, qc$kasp_id))
mi <- match(entries, mapping$phenotype_entry)
qi <- match(mapping$kasp_id[mi], qc$kasp_id)
stopifnot(!anyNA(mi), !anyNA(qi),
          all(mapping$sample_name[mi] == qc$sample_name[qi]))
identity_qc <- data.frame(entry = entries, kasp_id = mapping$kasp_id[mi],
  sample_name = mapping$sample_name[mi],
  kasp_missingness = qc$missing_call_proportion[qi],
  kasp_heterozygosity_valid_denominator = qc$heterozygous_call_proportion[qi])
stopifnot(!anyNA(identity_qc), all(identity_qc$kasp_missingness >= 0 & identity_qc$kasp_missingness <= 1),
          all(identity_qc$kasp_heterozygosity_valid_denominator >= 0 & identity_qc$kasp_heterozygosity_valid_denominator <= 1))
dir.create("outputs", showWarnings = FALSE)
dir.create("notes", showWarnings = FALSE)
membership <- list()
sets <- list()
ordered <- list()
for (s in scenarios) for (scale in scales) {
  z <- est[est$scenario == s & est$response_scale == scale, ]
  # Numeric ID is a display tie-breaker only; all cutoff ties are included.
  z <- z[order(-z$estimate, match(z$entry, entries)), ]
  z$display_order <- seq_len(nrow(z))
  ordered[[paste(s, scale)]] <- z
  for (k in c(5, 10)) {
    cutoff <- z$estimate[k]
    included <- z$estimate >= cutoff - tolerance
    key <- paste(s, scale, k, sep = "_")
    membership[[key]] <- data.frame(scenario = s, response_scale = scale,
      requested_size = k, entry = z$entry, estimate = z$estimate,
      display_order = z$display_order, cutoff_estimate = cutoff,
      absolute_tie_tolerance = tolerance, candidate = included,
      actual_set_size = sum(included))
    sets[[key]] <- z$entry[included]
  }
}
membership <- do.call(rbind, membership)
write.csv(do.call(rbind, ordered), "outputs/candidate_ordered_estimates.csv", row.names = FALSE)
write.csv(membership, "outputs/candidate_model_membership.csv", row.names = FALSE)

# Complete wide comparison: all 53 entries for BOTH requested set sizes.
# Means and intervals are labelled by response scale and scenario.
comparison <- list()
for (k in c(5, 10)) {
  x <- identity_qc
  x$requested_size <- k
  for (s in scenarios) for (scale in scales) {
    prefix <- paste(if (s == scenarios[1]) "A" else "B", scale, sep = "_")
    z <- est[est$scenario == s & est$response_scale == scale, ]
    z <- z[match(entries, z$entry), ]
    m <- membership[membership$scenario == s & membership$response_scale == scale & membership$requested_size == k, ]
    m <- m[match(entries, m$entry), ]
    for (field in c("estimate", "lower_95", "upper_95")) x[[paste(prefix, field, sep = "_")]] <- z[[field]]
    x[[paste0(prefix, "_candidate")]] <- m$candidate
    x[[paste0(prefix, "_actual_set_size")]] <- m$actual_set_size
  }
  candidate_columns <- grep("_candidate$", names(x), value = TRUE)
  x$models_including_entry <- rowSums(x[candidate_columns])
  comparison[[as.character(k)]] <- x
}
comparison <- do.call(rbind, comparison)
write.csv(comparison, "outputs/candidate_sensitivity_comparison.csv", row.names = FALSE)
write.csv(comparison[comparison$models_including_entry > 0, ],
          "outputs/candidate_sensitivity_candidate_details.csv", row.names = FALSE)

overlaps <- list()
for (k in c(5, 10)) {
  for (s in scenarios) {
    original <- sets[[paste(s, "original", k, sep = "_")]]
    sqrt_set <- sets[[paste(s, "sqrt", k, sep = "_")]]
    common <- intersect(original, sqrt_set)
    overlaps[[paste(s, k)]] <- data.frame(requested_size = k, scenario = s,
      original_size = length(original), sqrt_size = length(sqrt_set),
      overlap_count = length(common), union_count = length(union(original, sqrt_set)),
      overlap_ids = paste(common, collapse = ";"),
      original_only = paste(setdiff(original, sqrt_set), collapse = ";"),
      sqrt_only = paste(setdiff(sqrt_set, original), collapse = ";"))
  }
}
overlaps <- do.call(rbind, overlaps)
write.csv(overlaps, "outputs/candidate_scale_overlaps.csv", row.names = FALSE)

# Show only the union of candidates, arranged by numeric ID, not performance.
png("outputs/candidate_sensitivity_membership.png", width = 2000, height = 1500, res = 170)
par(mfrow = c(1, 2), mar = c(8, 5, 4, 1), oma = c(2, 0, 3, 0))
for (k in c(5, 10)) {
  x <- comparison[comparison$requested_size == k & comparison$models_including_entry > 0, ]
  cols <- c("A_original_candidate", "A_sqrt_candidate", "B_original_candidate", "B_sqrt_candidate")
  y <- rev(seq_len(nrow(x)))
  plot(NA, xlim = c(0.5, 4.5), ylim = c(0.5, nrow(x) + 0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = paste("Requested set size:", k))
  axis(1, at = 1:4, labels = c("A original", "A square root", "B original", "B square root"), las = 2)
  axis(2, at = y, labels = x$entry, las = 1)
  abline(h = y, col = "gray93")
  for (j in 1:4) {
    points(rep(j, nrow(x)), y, pch = 1, col = "gray80", cex = 1.5)
    points(rep(j, nrow(x))[x[[cols[j]]]], y[x[[cols[j]]]], pch = 16,
           col = if (j %in% c(1, 3)) "#286C87" else "#BD642A", cex = 1.5)
  }
}
mtext("Exploratory candidate membership — filled circles indicate inclusion", outer = TRUE, side = 3, line = 1, cex = 1.15)
mtext("Candidate unions in numeric ID order; membership does not establish superiority.", outer = TRUE, side = 1, cex = 0.9)
dev.off()
cat("Candidate sets (descending estimate, all cutoff ties included):\n")
set_lines <- character()
for (key in names(sets)) {
  line <- paste0(key, " (actual n=", length(sets[[key]]), "): ", paste(sets[[key]], collapse = ", "))
  cat(line, "\n")
  set_lines <- c(set_lines, paste0("- ", line))
}
cat("\nWithin-scenario scale overlaps:\n")
print(overlaps, row.names = FALSE)
all_four_lines <- vapply(c(5, 10), function(k) {
  x <- comparison[comparison$requested_size == k & comparison$models_including_entry == 4, ]
  paste0("Requested size ", k, ": ", nrow(x), " entries in all four models: ", paste(x$entry, collapse = ", "))
}, "")
cat(paste(all_four_lines, collapse = "\n"), "\n")
notes <- c(
  "# Exploratory candidate sensitivity", "",
  "Run `Rscript R/10_candidate_sensitivity.R` (base R). This is exploratory prioritization by grain weight per plant, not a declaration of statistically or genetically superior lines.",
  "The requested set sizes 5 and 10 are explicit exploratory choices, not biologically validated cutoffs. Within each of four models, entries are ordered by descending mean estimate. The kth estimate is the cutoff; include every estimate >= cutoff - 1e-10. This absolute tolerance is on that model's response scale and addresses numerical ties, not biological equivalence. Numeric entry ID only breaks display ties, never membership ties. Actual set sizes are reported.", "",
  "## Candidate sets", "", set_lines, "",
  "## Overlaps", "", capture.output(print(overlaps, row.names = FALSE)), "",
  all_four_lines, "",
  "## Estimates, intervals, and molecular QC", "",
  "The complete comparison CSV contains all 53 entries for each requested size, four membership indicators and actual sizes, and all four means and provisional pointwise intervals. Candidate details are also saved separately. Original-scale means target arithmetic mean grain weight equally averaged over blocks; square-root estimates target mean square-root grain weight. They are different quantities. Square-root estimates and limits are not squared or labelled arithmetic means. Units and sampling details remain unconfirmed.",
  "Candidate phenotype IDs are joined by explicit keys through sample_id_mapping.csv to kasp_sample_audit.csv. Unique IDs, full ID coverage, and agreement of sample names are checked. The phenotype-to-KASP correspondence remains user-confirmed, not independently field-verified. Every candidate's KASP missingness is missing calls divided by all 502 markers; heterozygosity is heterozygous calls divided by valid calls. The QC metrics come from the original unfiltered KASP audit. No candidate is excluded on molecular QC grounds.",
  "A uses all recorded observations; B omits only E-46/R-3 as sensitivity without declaring the observation erroneous. Agreement across models is sensitivity evidence, not independent validation and not a probability of superiority: all fits share the same trial and nearly all observations.",
  "Selecting candidates from the same data can exaggerate apparent performance (selection optimism). Intervals are pointwise and model-dependent, not simultaneous or selection-adjusted. Residual variance concerns remain; interval overlap is not a formal pairwise test. These results concern this trial, not broad adaptation, and further trials are needed.", "",
  "## Outputs", "",
  "- outputs/candidate_sensitivity_comparison.csv: complete wide comparison in numeric ID order within each requested size.",
  "- outputs/candidate_sensitivity_candidate_details.csv: union of candidates, with estimates, intervals, mapping, and molecular QC.",
  "- outputs/candidate_model_membership.csv: complete long membership report with cutoffs and tolerance.",
  "- outputs/candidate_ordered_estimates.csv: descending estimates within each model.",
  "- outputs/candidate_scale_overlaps.csv: original versus square-root overlap within each scenario.",
  "- outputs/candidate_sensitivity_membership.png: filled membership circles, in numeric entry order.",
  "No source data or previous results were changed. No marker associations, significance letters, or genetic-superiority labels are produced. No final set or model is selected.", "", R.version.string
)
writeLines(notes, "notes/10_candidate_sensitivity.md")
