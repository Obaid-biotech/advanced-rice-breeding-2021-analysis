# Run from the project folder: Rscript R/04_kasp_audit.R
# This audit uses readxl and base R, without changing or filtering source calls.
.libPaths(c(".R-library", .libPaths()))
library(readxl)
workbook <- "data/raw/NUYT-2021-KASP_GenotypeData.xlsx"
dir.create("outputs", showWarnings = FALSE)
dir.create("notes", showWarnings = FALSE)
save_report <- function(x, name) {
  write.csv(x, file.path("outputs", paste0("kasp_", name, ".csv")),
            row.names = FALSE, na = "NA")
}

# Read only the rectangle containing metadata and samples: columns A:SL.
# In Excel, column E is 5 and SL is 506. Summary rows/columns are excluded.
sheet <- as.matrix(read_excel(
  workbook, sheet = "Sheet1", range = "A1:SL55", col_names = FALSE,
  col_types = "text", trim_ws = FALSE, .name_repair = "minimal"
))
stopifnot(identical(dim(sheet), c(55L, 506L)))
sample_ids <- sheet[3:55, 1]
sample_names <- sheet[3:55, 2]
marker_names <- sheet[1, 5:506]
positions <- sheet[2, 5:506]
original_calls <- sheet[3:55, 5:506, drop = FALSE]
stopifnot(identical(dim(original_calls), c(53L, 502L)))
blank <- function(x) is.na(x) | trimws(x) == ""
duplicate_all <- function(x) {
  !blank(x) & (duplicated(x) | duplicated(x, fromLast = TRUE))
}

# Validate sample identity with an outer key join, never row order.
mapping <- read.csv("outputs/sample_id_mapping.csv", stringsAsFactors = FALSE)
stopifnot(all(c("kasp_id", "sample_name") %in% names(mapping)))
if (any(blank(mapping$kasp_id)) || anyDuplicated(mapping$kasp_id)) {
  stop("Mapping has missing or duplicate KASP IDs.")
}
samples <- data.frame(excel_row = 3:55, kasp_id = sample_ids,
                      sample_name = sample_names)
samples$missing_id <- blank(sample_ids)
samples$duplicate_id <- duplicate_all(sample_ids)
sample_check <- merge(samples, mapping[c("kasp_id", "sample_name")],
                      by = "kasp_id", all = TRUE, suffixes = c("_workbook", "_mapping"))
sample_check$unmatched <- is.na(sample_check$excel_row) |
  !sample_check$kasp_id %in% mapping$kasp_id
sample_check$name_mismatch <- is.na(sample_check$sample_name_workbook) |
  is.na(sample_check$sample_name_mapping) |
  sample_check$sample_name_workbook != sample_check$sample_name_mapping
save_report(sample_check, "sample_id_validation")
sample_ok <- !any(samples$missing_id | samples$duplicate_id) &&
  !any(sample_check$unmatched | sample_check$name_mismatch) && nrow(mapping) == 53

# Metadata problems are reported without removing markers. Position syntax
# expected here is chromosome:coordinate; syntax checks do not verify a genome build.
metadata <- data.frame(
  excel_column = 5:506, marker_name = marker_names, position = positions,
  missing_marker_name = blank(marker_names),
  duplicate_marker_name = duplicate_all(marker_names),
  missing_position = blank(positions),
  unexpected_position_format = !blank(positions) &
    !grepl("^([1-9]|1[0-2]):[1-9][0-9]*$", positions),
  repeated_position = duplicate_all(positions)
)
save_report(metadata, "marker_metadata_audit")
save_report(metadata[metadata$repeated_position, ], "repeated_positions")

# Preserve original text. Empty/whitespace-only cells and the three explicit
# missing labels are missing; other text must be an exact A/C/G/T allele pair.
missing <- is.na(original_calls) | trimws(original_calls) == "" |
  original_calls %in% c("?", "Uncallable", "Bad")
dim(missing) <- dim(original_calls)
valid <- matrix(grepl("^[ACGT]:[ACGT]$", original_calls), nrow = 53)
unexpected <- !missing & !valid
labels <- as.vector(original_calls)
labels[is.na(labels)] <- "<blank cell>"
labels[!is.na(labels) & trimws(labels) == ""] <- "<blank cell>"
frequencies <- as.data.frame(table(labels), stringsAsFactors = FALSE)
names(frequencies) <- c("original_call_string", "frequency")
save_report(frequencies, "original_call_frequencies")

# Record unexpected cells by Excel coordinates and keys for inspection.
where <- which(unexpected, arr.ind = TRUE)
unexpected_records <- data.frame(
  excel_row = where[, 1] + 2L, excel_column = where[, 2] + 4L,
  kasp_id = sample_ids[where[, 1]], marker_name = marker_names[where[, 2]],
  original_call = original_calls[where]
)
save_report(unexpected_records, "unexpected_calls")

# Sort each valid pair alphabetically: G:A and A:G both become A:G.
# Unexpected calls are unscorable and become NA only in this derived copy.
analysis_calls <- matrix(NA_character_, nrow = 53, ncol = 502)
pairs <- strsplit(original_calls[valid], ":", fixed = TRUE)
analysis_calls[valid] <- vapply(pairs, function(x) paste(sort(x), collapse = ":"), "")
heterozygous <- valid & matrix(
  substr(original_calls, 1, 1) != substr(original_calls, 3, 3), nrow = 53
)
heterozygous[is.na(heterozygous)] <- FALSE
valid_per_sample <- rowSums(valid)
sample_audit <- data.frame(
  kasp_id = sample_ids, sample_name = sample_names,
  total_markers = ncol(original_calls), missing_calls = rowSums(missing),
  unexpected_calls = rowSums(unexpected), valid_calls = valid_per_sample,
  missing_call_proportion = rowMeans(missing),
  heterozygous_calls = rowSums(heterozygous),
  heterozygous_call_proportion = ifelse(valid_per_sample > 0,
    rowSums(heterozygous) / valid_per_sample, NA_real_)
)
alleles <- lapply(seq_len(ncol(original_calls)), function(j) {
  sort(unique(unlist(strsplit(original_calls[valid[, j], j], ":", fixed = TRUE))))
})
allele_counts <- lengths(alleles)
marker_audit <- data.frame(
  excel_column = 5:506, marker_name = marker_names, position = positions,
  total_samples = nrow(original_calls), missing_calls = colSums(missing),
  unexpected_calls = colSums(unexpected), valid_calls = colSums(valid),
  missing_call_proportion = colMeans(missing),
  observed_alleles = vapply(alleles, paste, "", collapse = "/"),
  number_observed_alleles = allele_counts,
  classification = ifelse(allele_counts == 0, "no_valid_calls",
    ifelse(allele_counts == 1, "monomorphic",
      ifelse(allele_counts == 2, "biallelic", "multiallelic"))),
  entirely_missing = colSums(missing) == nrow(original_calls)
)
save_report(sample_audit, "sample_audit")
save_report(marker_audit, "marker_audit")

# Long copies keep every cell, the original call, and its derived counterpart.
# This joins no phenotype measurements and retains every sample and marker.
call_copies <- data.frame(
  kasp_id = rep(sample_ids, times = 502),
  excel_column = rep(5:506, each = 53),
  marker_name = rep(marker_names, each = 53),
  original_call = as.vector(original_calls),
  analysis_call = as.vector(analysis_calls),
  call_status = ifelse(as.vector(missing), "missing",
    ifelse(as.vector(valid), "valid", "unexpected"))
)
save_report(call_copies, "original_and_analysis_calls")

summary <- data.frame(
  measure = c("samples", "markers", "total_cells", "missing_calls", "unexpected_calls",
    "monomorphic_markers", "biallelic_markers", "multiallelic_markers",
    "entirely_missing_markers", "markers_without_valid_calls",
    "missing_marker_names", "markers_with_duplicate_names", "missing_positions",
    "unexpected_position_formats", "markers_with_repeated_positions",
    "distinct_repeated_positions"),
  count = c(53, 502, length(original_calls), sum(missing), sum(unexpected),
    sum(allele_counts == 1), sum(allele_counts == 2), sum(allele_counts > 2),
    sum(marker_audit$entirely_missing), sum(allele_counts == 0),
    sum(metadata$missing_marker_name), sum(metadata$duplicate_marker_name),
    sum(metadata$missing_position), sum(metadata$unexpected_position_format),
    sum(metadata$repeated_position), length(unique(positions[metadata$repeated_position])))
)
save_report(summary, "audit_summary")
versions <- data.frame(component = c("R", "readxl"),
                       version = c(R.version.string, as.character(packageVersion("readxl"))))
save_report(versions, "software_versions")
cat("Sample identity validation passed:", sample_ok, "\n")
print(summary, row.names = FALSE)
cat("\nEvery distinct original call string and its frequency:\n")
print(frequencies, row.names = FALSE)
cat("\nSample missing-call proportion range:", range(sample_audit$missing_call_proportion),
    "\nMarker missing-call proportion range:", range(marker_audit$missing_call_proportion),
    "\nSample heterozygous proportion range:",
    range(sample_audit$heterozygous_call_proportion, na.rm = TRUE), "\n")
print(versions, row.names = FALSE)

# Generate notes from this run, keeping recorded results reproducible.
notes <- c(
  "# 2021 KASP audit", "",
  "Run `Rscript R/04_kasp_audit.R` from the project folder.", "",
  "Sheet1 rows 3–55 and columns 5–506 (E:SL) contain 53 samples and 502 markers.",
  "Marker names come from row 1 and positions from row 2. Headers, sample metadata, and summary cells are excluded from the call matrix.",
  paste("Sample IDs and Sample Name text agree with sample_id_mapping.csv:", sample_ok),
  "The identity correspondence remains user-confirmed, not independently verified from field records.", "",
  "## Results", "",
  paste(summary$measure, summary$count, sep = ": "), "",
  paste("Sample missing proportion range:", paste(range(sample_audit$missing_call_proportion), collapse = " to ")),
  paste("Marker missing proportion range:", paste(range(marker_audit$missing_call_proportion), collapse = " to ")),
  paste("Sample heterozygous proportion range:", paste(range(sample_audit$heterozygous_call_proportion, na.rm = TRUE), collapse = " to ")), "",
  "Every original call string and its frequency:", "",
  paste(frequencies$original_call_string, frequencies$frequency, sep = ": "), "",
  "## Definitions and denominators", "",
  "Missing calls are ?, Uncallable, Bad, and blank cells (including whitespace-only cells). Sample missing proportion divides missing cells by all 502 markers. Marker missing proportion divides missing cells by all 53 samples.",
  "Valid calls are exact A/C/G/T pairs separated by a colon. Unexpected strings are reported separately, not silently counted as the defined missing labels; they are unscorable in the derived copy.",
  "Sample heterozygous proportion divides valid pairs with two different alleles by all valid pairs for that sample. Missing and unexpected cells are excluded from this denominator. If there are no valid calls, the proportion is NA.",
  "Observed allele counts use distinct A/C/G/T letters across valid pairs only. One allele means monomorphic, two biallelic, and more than two multiallelic. Entirely missing means all 53 cells have a defined missing label; markers with no valid calls are also counted separately.",
  "The derived analysis copy sorts valid pairs alphabetically, so G:A and A:G both become A:G. Original text is retained alongside it. Repeated positions are reported without removing markers. Position syntax checks do not verify a reference genome or genome build.", "",
  "## Outputs", "",
  "All audit reports have the outputs/kasp_ prefix: sample_id_validation, marker_metadata_audit, repeated_positions, original_call_frequencies, unexpected_calls, sample_audit, marker_audit, original_and_analysis_calls, audit_summary, and software_versions (all CSV). Empty reports retain column headers.",
  "The original_and_analysis_calls report retains every cell and its status. Original blank cells are written as NA; their original cell locations are retained through sample IDs and Excel marker columns.",
  "No samples or markers were filtered, calls imputed, marker calls joined to phenotype measurements, or marker–trait associations tested.", "",
  paste(versions$component, versions$version, sep = ": ")
)
writeLines(notes, "notes/04_kasp_audit.md")
