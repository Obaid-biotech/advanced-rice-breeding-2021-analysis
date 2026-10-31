# Run from the project folder: Rscript R/03_sample_mapping.R
# Only readxl and base R are used; no marker calls are imported or merged.
.libPaths(c(".R-library", .libPaths()))
library(readxl)

phenotype <- read.csv("outputs/phenotype_2021_imported.csv",
                      stringsAsFactors = FALSE, check.names = FALSE)
if (!"entry" %in% names(phenotype)) stop("Phenotype CSV lacks entry column.")

# Replications repeat entry IDs intentionally. Reduce to distinct entries,
# but first count missing labels in the original records.
missing_id <- function(x) is.na(x) | trimws(x) == ""
cat("Missing phenotype entry labels:", sum(missing_id(phenotype$entry)), "\n")
entries <- unique(phenotype$entry)
cat("Phenotype records:", nrow(phenotype), "; distinct entries:", length(entries), "\n")
cat("Repeated phenotype labels across replications (expected):",
    sum(duplicated(phenotype$entry)), "\n")

# Inspection of Sheet1 shows actual samples in rows 3:55, with NIBG IDs
# in column A and Sample Name in column B. This explicit sample range
# excludes the header, Indica position row, and later summary rows.
kasp <- as.data.frame(read_excel(
  "data/raw/Advanced_lines-2021-KASP_GenotypeData.xlsx", sheet = "Sheet1",
  range = "A3:B55", col_names = c("kasp_id", "sample_name"),
  col_types = "text", trim_ws = FALSE
))
cat("KASP sample rows:", nrow(kasp), "\n")
cat("Missing KASP IDs:", sum(missing_id(kasp$kasp_id)), "\n")
cat("Missing sample names:", sum(missing_id(kasp$sample_name)), "\n")
cat("Duplicate KASP IDs:", sum(duplicated(kasp$kasp_id)), "\n")
cat("Duplicate sample names:", sum(duplicated(kasp$sample_name)), "\n")
if (any(missing_id(entries)) || any(missing_id(kasp$kasp_id)) ||
    any(missing_id(kasp$sample_name)) || anyDuplicated(kasp$kasp_id) ||
    anyDuplicated(kasp$sample_name)) {
  stop("Missing or duplicate identifiers prevent a one-to-one mapping.")
}

# The USER confirmed E-n corresponds to NYUT_21_n for n = 1,...,53.
# This is not independent verification from field records.
# Validate exact formats so unexpected spelling is never silently changed.
if (!all(grepl("^E-[1-9][0-9]*$", entries)) ||
    !all(grepl("^NYUT_21_[1-9][0-9]*$", kasp$kasp_id))) {
  stop("Unexpected ID format; exact NYUT_21_ spelling is required.")
}
phenotype_suffix <- as.integer(sub("^E-", "", entries))
kasp_suffix <- as.integer(sub("^NYUT_21_", "", kasp$kasp_id))
if (anyDuplicated(phenotype_suffix) || anyDuplicated(kasp_suffix)) {
  stop("Duplicate numeric suffixes prevent one-to-one correspondence.")
}
cat("Missing expected phenotype IDs:\n")
print(setdiff(paste0("E-", 1:53), entries))
cat("Missing expected KASP IDs:\n")
print(setdiff(paste0("NYUT_21_", 1:53), kasp$kasp_id))

# match() searches for the numeric key, never for the row position.
index <- match(phenotype_suffix, kasp_suffix)
unmatched_entries <- entries[is.na(index)]
unmatched_kasp <- kasp$kasp_id[!kasp_suffix %in% phenotype_suffix]
cat("Unmatched phenotype entries:\n")
print(unmatched_entries)
cat("Unmatched KASP IDs:\n")
print(unmatched_kasp)
one_to_one <- length(entries) == 53 && nrow(kasp) == 53 &&
  setequal(phenotype_suffix, 1:53) && setequal(kasp_suffix, 1:53) &&
  length(unmatched_entries) == 0 && length(unmatched_kasp) == 0 &&
  !anyDuplicated(index)
cat("One-to-one correspondence for all 53 entries:", one_to_one, "\n")
if (!one_to_one) stop("Mapping validation failed; no mapping saved.")

mapping <- data.frame(
  phenotype_entry = entries, kasp_id = kasp$kasp_id[index],
  sample_name = kasp$sample_name[index], stringsAsFactors = FALSE
)
# Sorting is only for readable output, after matching by key.
mapping <- mapping[order(phenotype_suffix), ]
dir.create("outputs", showWarnings = FALSE)
write.csv(mapping, "outputs/sample_id_mapping.csv", row.names = FALSE)
cat("\nFirst six mapping rows:\n")
print(head(mapping), row.names = FALSE)
cat("Saved outputs/sample_id_mapping.csv\n")
