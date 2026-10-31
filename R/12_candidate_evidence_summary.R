# Run from project folder: Rscript R/12_candidate_evidence_summary.R
# Base R only. Existing scripts 03, 04, 10, 11 and CSV structures were inspected.
read_output <- function(f) read.csv(file.path("outputs", f), stringsAsFactors = FALSE)
map <- read_output("sample_id_mapping.csv")
qc <- read_output("candidate_marker_profiles_sample_qc.csv")
candidates <- read_output("candidate_sensitivity_candidate_details.csv")
candidates <- candidates[candidates$requested_size == 10, ]
five <- read_output("candidate_model_membership.csv")
five <- five[five$requested_size == 5, ]
phen <- read_output("phenotype_2021_imported.csv")
calls <- read_output("kasp_original_and_analysis_calls.csv")
markers <- read_output("kasp_marker_audit.csv")
factor_flags <- read_output("yield_factor_discrepancies.csv")
sheet_flags <- read_output("yield_between_sheet_discrepancies.csv")
entries <- paste0("E-", 1:53)
candidates <- candidates[order(match(candidates$entry, entries)), ]
stopifnot(nrow(candidates) == 10, !anyDuplicated(candidates$entry),
          !anyDuplicated(map$phenotype_entry), !anyDuplicated(map$kasp_id),
          !anyDuplicated(qc$phenotype_entry), nrow(markers) == 502,
          !anyDuplicated(markers$excel_column),
          !anyDuplicated(calls[c("kasp_id", "excel_column")]),
          !anyDuplicated(phen[c("entry", "replication")]))
mi <- match(candidates$entry, map$phenotype_entry)
qi <- match(candidates$entry, qc$phenotype_entry)
stopifnot(!anyNA(mi), !anyNA(qi),
          all(map$kasp_id[mi] == candidates$kasp_id),
          all(map$kasp_id[mi] == qc$kasp_id[qi]),
          all(map$sample_name[mi] == qc$sample_name[qi]),
          all(map$kasp_id[mi] == paste0("Advacned_21_", sub("E-", "", candidates$entry))))

# Minimal read-only extraction of A:D metadata from this inspected XLSX layout.
# XLSX is a ZIP of XML files. This handles its shared-string cells, not arbitrary
# Excel files. No external R package or workbook modification is needed.
workbook <- "data/raw/Advanced-2021-KASP_GenotypeData.xlsx"
read_xml <- function(member) {
  con <- unz(workbook, member); on.exit(close(con))
  paste(readLines(con, warn = FALSE), collapse = "")
}
extract <- function(pattern, x) regmatches(x, gregexpr(pattern, x, perl = TRUE))[[1]]
decode <- function(x) {
  for (entity in c("&lt;", "&gt;", "&quot;", "&apos;", "&amp;")) {
    replacement <- c("<", ">", '"', "'", "&")[match(entity, c("&lt;", "&gt;", "&quot;", "&apos;", "&amp;"))]
    x <- gsub(entity, replacement, x, fixed = TRUE)
  }
  x
}
strings <- vapply(extract("<si(?:\\s[^>]*)?>.*?</si>", read_xml("xl/sharedStrings.xml")),
  function(x) decode(paste(gsub("<[^>]+>", "", extract("<t(?:\\s[^>]*)?>.*?</t>", x)), collapse = "")), "")
sheet <- read_xml("xl/worksheets/sheet1.xml")
cell <- function(address) {
  x <- extract(paste0('<c\\b[^>]*\\br="', address, '"[^>]*>.*?</c>'), sheet)
  if (length(x) != 1) stop("Missing or ambiguous metadata cell: ", address)
  if (!grepl('t="s"', x, fixed = TRUE)) stop("Unexpected metadata cell type.")
  value <- extract("<v>[^<]+</v>", x)
  stopifnot(length(value) == 1)
  strings[as.integer(gsub("<[^>]+>", "", value)) + 1L]
}
stopifnot(cell("A1") == "NIBG", cell("B1") == "Sample Name",
          cell("C1") == "Master plate", cell("D1") == "Master well")
metadata <- as.data.frame(t(vapply(3:55, function(r) vapply(LETTERS[1:4],
  function(col) cell(paste0(col, r)), ""), character(4))), stringsAsFactors = FALSE)
names(metadata) <- c("kasp_id", "sample_name", "master_plate", "master_well")
stopifnot(!anyDuplicated(metadata$kasp_id), setequal(metadata$kasp_id, map$kasp_id),
          all(metadata$sample_name[match(map$kasp_id, metadata$kasp_id)] == map$sample_name))
ids <- map$kasp_id[match(c("E-42", "E-45"), map$phenotype_entry)]
meta_pair <- metadata[match(ids, metadata$kasp_id), ]
meta_pair$phenotype_entry <- c("E-42", "E-45")
stopifnot(length(unique(meta_pair$sample_name)) == 2)
write.csv(meta_pair, "outputs/e42_e45_laboratory_metadata.csv", row.names = FALSE)

pair_calls <- lapply(ids, function(id) {
  z <- calls[calls$kasp_id == id, ]; z <- z[match(markers$excel_column, z$excel_column), ]
  stopifnot(nrow(z) == 502, !anyNA(z$excel_column))
  original <- z$original_call
  valid <- grepl("^[ACGT]:[ACGT]$", original)
  missing <- is.na(original) | trimws(original) == "" | original %in% c("?", "Bad", "Uncallable")
  stopifnot(all(valid | missing))
  canonical <- rep(NA_character_, 502)
  canonical[valid] <- vapply(strsplit(original[valid], ":", fixed = TRUE),
                             function(x) paste(sort(x), collapse = ":"), "")
  stopifnot(identical(is.na(canonical), is.na(z$analysis_call)),
            all(canonical[valid] == z$analysis_call[valid]))
  z
})
a <- pair_calls[[1]]; b <- pair_calls[[2]]
va <- !is.na(a$analysis_call); vb <- !is.na(b$analysis_call)
category <- rep(NA_character_, 502)
category[va & vb] <- ifelse(a$analysis_call[va & vb] == b$analysis_call[va & vb], "matching_valid", "mismatching_valid")
category[!va & vb] <- "missing_E42_only"
category[va & !vb] <- "missing_E45_only"
category[!va & !vb] <- "missing_both"
levels_category <- c("matching_valid", "mismatching_valid", "missing_E42_only", "missing_E45_only", "missing_both")
counts <- data.frame(category = levels_category,
                     count = as.integer(table(factor(category, levels = levels_category))))
stopifnot(!anyNA(category), sum(counts$count) == 502,
          sum(va & vb) == sum(counts$count[1:2]))
marker_comparison <- data.frame(excel_column = markers$excel_column,
  marker_name = markers$marker_name, position = markers$position,
  E42_original_call = a$original_call, E45_original_call = b$original_call,
  E42_canonical_call = a$analysis_call, E45_canonical_call = b$analysis_call,
  category = category)
write.csv(counts, "outputs/e42_e45_marker_category_counts.csv", row.names = FALSE)
write.csv(data.frame(valid_in_both = sum(va & vb), category_total = sum(counts$count), expected_markers = 502),
          "outputs/e42_e45_marker_total_validation.csv", row.names = FALSE)
write.csv(marker_comparison, "outputs/e42_e45_all_marker_comparison.csv", row.names = FALSE)
write.csv(marker_comparison[category == "mismatching_valid", ], "outputs/e42_e45_mismatching_markers.csv", row.names = FALSE)

# Original phenotype values, explicitly side by side within each replication.
traits <- c("plant_height", "tillers", "panicle_length", "filled_grains",
            "unfilled_grains", "grain_weight_100", "plant_weight", "recorded_yield")
stopifnot(all(vapply(phen[traits], is.numeric, TRUE)))
side <- data.frame(replication = c("R-1", "R-2", "R-3"))
for (entry in c("E-42", "E-45")) {
  z <- phen[phen$entry == entry, ]
  z <- z[match(side$replication, z$replication), ]
  stopifnot(nrow(z) == 3, !anyNA(z$entry))
  for (trait in traits) side[[paste(entry, trait, sep = "_")]] <- z[[trait]]
}
write.csv(side, "outputs/e42_e45_phenotypes_side_by_side.csv", row.names = FALSE)

# Retain both scenarios and scales from script 09, with explicit column labels.
evidence <- candidates[c("entry", "kasp_id", "sample_name",
  grep("^(A|B)_(original|sqrt)_(estimate|lower_95|upper_95)$", names(candidates), value = TRUE))]
for (model in unique(paste(five$scenario, five$response_scale))) {
  z <- five[paste(five$scenario, five$response_scale) == model, ]
  index <- match(evidence$entry, z$entry)
  stopifnot(!anyNA(index))
  evidence[[paste0("five_entry_set_", gsub(" ", "_", model))]] <- z$candidate[index]
}
other_traits <- setdiff(traits, c("plant_weight", "recorded_yield"))
for (trait in other_traits) {
  evidence[[paste0("descriptive_mean_", trait)]] <- vapply(evidence$entry,
    function(e) mean(phen[[trait]][phen$entry == e], na.rm = TRUE), 0)
}
# recorded_yield is deliberately separate, with its unresolved interpretation.
evidence$descriptive_mean_recorded_yield_unresolved <- vapply(evidence$entry,
  function(e) mean(phen$recorded_yield[phen$entry == e]), 0)
for (field in setdiff(names(qc), c("phenotype_entry", "kasp_id", "sample_name"))) evidence[[field]] <- qc[[field]][qi]
evidence$factor_discrepancy_records <- vapply(evidence$entry, function(e) sum(factor_flags$entry == e), 0L)
evidence$between_sheet_yield_discrepancy_records <- vapply(evidence$entry, function(e) sum(sheet_flags$entry == e), 0L)
evidence$observed_zero_measurements <- vapply(evidence$entry,
  function(e) sum(phen[phen$entry == e, traits] == 0, na.rm = TRUE), 0L)
evidence$identity_review_E42_E45 <- evidence$entry %in% c("E-42", "E-45")
evidence$measurement_details <- "plant_weight = grain weight per plant (user-confirmed); unit and sampling unconfirmed; recorded_yield conversion unresolved"
evidence$interval_limitations <- "pointwise model-dependent mean CIs; residual variance concerns; not simultaneous or selection-adjusted"
evidence$evidence_supports <- "membership in exploratory ten-entry sets across four fits of this 2021 trial; descriptive traits and observed targeted-panel genotype information"
evidence$remaining_uncertainties <- "selection optimism; sampling and units; model assumptions; broader adaptation and breeding objective; no genetic superiority established"
write.csv(evidence, "outputs/candidate_evidence_summary.csv", row.names = FALSE)

candidate_notes <- vapply(seq_len(nrow(evidence)), function(i) {
  x <- evidence[i, ]
  extra <- character()
  if (x$factor_discrepancy_records > 0) extra <- c(extra, "recorded_yield factor discrepancy requires source review")
  if (x$between_sheet_yield_discrepancy_records > 0) extra <- c(extra, "Combined versus replication recorded_yield discrepancy requires source review")
  if (x$identity_review_E42_E45) extra <- c(extra, "E-42/E-45 observed-call agreement requires identity provenance review, without an identity conclusion")
  if (x$exceeds_exploratory_sample_missingness_10pct) extra <- c(extra, "exceeds exploratory 10% sample missingness limit (retained)")
  paste0("- ", x$entry, ": original A grain-weight estimate ", sprintf("%.3f", x$A_original_estimate),
    " (pointwise 95% CI ", sprintf("%.3f to %.3f", x$A_original_lower_95, x$A_original_upper_95),
    "); missingness ", sprintf("%.2f%%", 100*x$missing_call_proportion),
    ", heterozygosity among valid calls ", sprintf("%.2f%%",100*x$heterozygous_call_proportion),
    ". Supports exploratory follow-up within this trial only. ",
    if (length(extra)) paste0("Specific cautions: ", paste(extra, collapse = "; "), ". ") else "",
    "Units, sampling, model assumptions, selection optimism, and wider adaptation remain uncertain; no trait-weighting objective is defined.")
}, "")
notes <- c("# Candidate evidence summary", "",
  "Run `Rscript R/12_candidate_evidence_summary.R` (base R only). Existing scripts and output structures were inspected. Prior data and outputs are preserved. Validated key joins use the user-confirmed phenotype E-n to Advanced_21_n correspondence; this is not independent field-record verification.", "",
  "## E-42 and E-45: all 502 original markers", "",
  paste("Valid calls in both:", sum(va & vb)),
  capture.output(print(counts, row.names = FALSE)),
  "Matching and mismatching valid counts partition valid-in-both; do not add valid-in-both a second time. The five disjoint categories total exactly 502. Missing calls remain ?, Bad, Uncallable, or blank; canonical reversed pairs agree, while heterozygous and homozygous genotypes differ. All mismatches and their original calls are exported, even if the resulting CSV contains headers only.", "",
  capture.output(print(meta_pair, row.names = FALSE)),
  "Laboratory metadata are extracted read-only from Sheet1 A:D via the workbook ZIP/XML using base R and checked against the mapping. Different sample names or wells document different labels/locations, not independently established biological identity.", "",
  "Original phenotype values side by side (plant_weight means grain weight per plant, unit unconfirmed):", "",
  capture.output(print(side, row.names = FALSE)),
  "recorded_yield is a separate unresolved measurement. The 135.85 relationship is unverified; audit 02 discrepancies are carried into the evidence table and are not corrected.", "",
  "## Possible explanations and needed evidence", "",
  "Possible explanations include distinct materials not distinguished by this targeted marker set, shared ancestry or closely related material, a seed/sample handling or labelling issue, or a genotyping/data transfer issue. These are alternatives to investigate, not conclusions. Neither identical observed panel calls nor different phenotype measurements establishes biological identity or its absence.",
  "Review field entry books and randomization maps, seed inventories and pedigree/accession records, collection labels and chain of custody, DNA extraction logs, plate maps, assay cluster/raw fluorescence files, and export/transcription history. Independently collect traceable plants or seed lots, repeat DNA extraction and genotyping with controls, and consider additional informative genome-distributed markers or sequencing. Replicate sampling within entries may help distinguish within-entry variation from handling problems. No identities are corrected and no crosses are recommended.", "",
  "## Ten-candidate evidence", "", candidate_notes, "",
  "candidate_evidence_summary.csv contains membership in each model's five-entry set; original-scale estimates and provisional mean-response pointwise CIs for A and B; square-root estimates and CIs on their own scale; descriptive means of the other six traits; a separately labelled unresolved recorded_yield mean; original molecular QC metrics and exploratory flags; and specific audit counts.",
  "Other trait means are descriptions, not favorable directions. No direction, weight, or composite score is assigned without a defined breeding objective. All ten remain retained, including E-37. Molecular missingness uses all 502 markers; heterozygosity uses valid calls. Selection on these same trial data may exaggerate performance. Agreement across the four related fits is sensitivity evidence, not independent validation or probability of superiority.",
  "A includes all observations, B omits only E-46/R-3 as previously defined sensitivity. Square-root means are not squared into arithmetic means. Confidence intervals are pointwise, model-dependent, not simultaneous, not prediction intervals, and not adjusted for candidate selection. Residual variance concerns persist. Diagnostics cannot verify randomization or independence; units and sampling details remain unresolved. These results concern this trial, not broad adaptation, and further trials are needed.", "",
  "## Outputs", "",
  "outputs/candidate_evidence_summary.csv; e42_e45_marker_category_counts.csv; e42_e45_marker_total_validation.csv; e42_e45_all_marker_comparison.csv; e42_e45_mismatching_markers.csv; e42_e45_laboratory_metadata.csv; e42_e45_phenotypes_side_by_side.csv.",
  "No candidates removed, identities corrected, marker associations tested, significance letters assigned, or crosses recommended.", "", R.version.string)
writeLines(notes, "notes/12_candidate_evidence_summary.md")
cat("All-marker category counts (total", sum(counts$count), "); valid in both:", sum(va & vb), "\n")
print(counts, row.names = FALSE)
print(meta_pair, row.names = FALSE)
print(side, row.names = FALSE)
cat("Evidence table:", nrow(evidence), "candidates; key and category checks passed.\n")
