# Run from the project folder: Rscript R/02_yield_consistency.R
# Use the existing project library. Only readxl and base R are needed.
.libPaths(c(".R-library", .libPaths()))
library(readxl)
tolerance <- 0.000001
factor_value <- 135.85

phenotype <- read.csv("outputs/phenotype_2021_imported.csv",
                      stringsAsFactors = FALSE, check.names = FALSE)
expected_names <- c(
  "entry", "plant_height", "tillers", "panicle_length", "filled_grains",
  "unfilled_grains", "grain_weight_100", "plant_weight", "recorded_yield",
  "entry_original", "replication"
)
cat("Imported CSV column names:\n")
print(names(phenotype))
if (!identical(names(phenotype), expected_names)) {
  stop("CSV column names or their order differ from the phenotype import.")
}
stopifnot(is.numeric(phenotype$plant_weight),
          is.numeric(phenotype$recorded_yield))

# This factor is an unverified numerical relationship, NOT a validated
# conversion to yield per hectare. Calculations never replace observations.
phenotype$expected_from_factor <- factor_value * phenotype$plant_weight
phenotype$difference <- phenotype$recorded_yield - phenotype$expected_from_factor
report_columns <- c("entry", "replication", "plant_weight", "recorded_yield",
                    "expected_from_factor", "difference")
# abs() ignores the sign; > tests whether the magnitude exceeds tolerance.
factor_flag <- !is.na(phenotype$difference) &
  abs(phenotype$difference) > tolerance
factor_discrepancies <- phenotype[factor_flag, report_columns]
cat("\nFactor discrepancies:", nrow(factor_discrepancies), "\n")
print(factor_discrepancies, row.names = FALSE, digits = 15)
cat("Rows with an uncheckable factor difference:",
    sum(is.na(phenotype$difference)), "\n")

# == tests equality, | means OR, and & means AND.
# A logical vector inside [rows, columns] selects rows where it is TRUE.
zero_flag <- (!is.na(phenotype$plant_weight) & phenotype$plant_weight == 0) |
  (!is.na(phenotype$recorded_yield) & phenotype$recorded_yield == 0)
zero_records <- phenotype[zero_flag, report_columns]
cat("\nRecords with zero plant weight or recorded yield:\n")
print(zero_records, row.names = FALSE, digits = 15)
cat("A zero-zero pair satisfies 0 = 135.85 * 0, but the equation cannot\n",
    "tell whether zeros represent biological observations or recording issues.\n",
    "Their meaning requires clarification from field/source records.\n", sep = "")

combined <- as.data.frame(read_excel(
  "data/raw/NUYT-2021 Orignal .xlsx", sheet = "Combined", range = "A3:E55",
  col_names = c("entry", "yield_r1", "yield_r2", "yield_r3", "mean"),
  col_types = c("text", rep("numeric", 4)), trim_ws = FALSE
))
combined$entry_original <- combined$entry
combined$entry <- trimws(combined$entry)

# Validate IDs before matching: duplicate keys would make comparisons ambiguous.
expected_ids <- paste0("E-", 1:53)
cat("\nCombined rows:", nrow(combined),
    "; unique nonmissing IDs:", length(unique(na.omit(combined$entry))), "\n")
cat("Combined missing expected IDs:\n")
print(setdiff(expected_ids, combined$entry))
cat("Combined unexpected IDs:\n")
print(setdiff(combined$entry, expected_ids))
bad_id <- is.na(combined$entry) | combined$entry == "" |
  duplicated(combined$entry) | duplicated(combined$entry, fromLast = TRUE)
cat("Combined missing, blank, or duplicate ID records:\n")
print(combined[bad_id, ], row.names = FALSE)
if (any(bad_id)) stop("Resolve ambiguous Combined IDs before key comparison.")
if (anyNA(phenotype[c("entry", "replication")]) ||
    any(duplicated(phenotype[c("entry", "replication")])) ||
    any(phenotype$entry == "" | phenotype$replication == "")) {
  stop("Phenotype entry/replication keys must be present and unique.")
}

# Convert the three Combined yield columns to one row per entry/replication.
combined_long <- rbind(
  data.frame(entry = combined$entry, replication = "R-1", combined_yield = combined$yield_r1),
  data.frame(entry = combined$entry, replication = "R-2", combined_yield = combined$yield_r2),
  data.frame(entry = combined$entry, replication = "R-3", combined_yield = combined$yield_r3)
)
combined_long$in_combined <- TRUE
replication_yields <- phenotype[c("entry", "replication", "recorded_yield")]
replication_yields$in_replication <- TRUE
# all = TRUE retains unmatched keys from BOTH sources. merge() joins by key,
# so worksheet row order has no effect on the comparison.
comparison <- merge(combined_long, replication_yields,
                    by = c("entry", "replication"), all = TRUE)
comparison$difference <- comparison$combined_yield - comparison$recorded_yield
comparison$status <- "matched"
comparison$status[is.na(comparison$in_combined)] <- "missing_from_combined"
comparison$status[is.na(comparison$in_replication)] <- "missing_from_replication"
matched <- comparison$status == "matched"
comparison$status[matched & is.na(comparison$difference)] <- "missing_yield"
comparison$status[matched & !is.na(comparison$difference) &
                      abs(comparison$difference) > tolerance] <- "yield_difference"
unmatched <- !matched
cat("\nUnmatched entry/replication keys:", sum(unmatched), "\n")
print(comparison[unmatched, ], row.names = FALSE)
sheet_discrepancies <- comparison[comparison$status != "matched", ]
cat("\nBetween-sheet discrepancies (Combined minus replication):\n")
print(sheet_discrepancies, row.names = FALSE, digits = 15)

# Do not remove missing values: a missing input makes the mean uncheckable.
combined$expected_mean <- rowMeans(combined[c("yield_r1", "yield_r2", "yield_r3")])
combined$mean_difference <- combined$mean - combined$expected_mean
mean_flag <- is.na(combined$mean_difference) |
  abs(combined$mean_difference) > tolerance
mean_discrepancies <- combined[mean_flag, ]
cat("\nCombined mean discrepancies or uncheckable means:",
    nrow(mean_discrepancies), "\n")
print(mean_discrepancies, row.names = FALSE, digits = 15)

# Empty reports still have headers. These are new audit outputs, not corrections.
dir.create("outputs", showWarnings = FALSE)
write.csv(factor_discrepancies, "outputs/yield_factor_discrepancies.csv", row.names = FALSE)
write.csv(zero_records, "outputs/yield_zero_records.csv", row.names = FALSE)
write.csv(sheet_discrepancies, "outputs/yield_between_sheet_discrepancies.csv", row.names = FALSE)
write.csv(mean_discrepancies, "outputs/yield_combined_mean_discrepancies.csv", row.names = FALSE)
