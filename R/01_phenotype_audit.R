# Audit the three 2021 replications without correcting measurements.
# Run from the project folder: Rscript R/01_phenotype_audit.R

# Look for readxl in the existing project package library first.
.libPaths(c(".R-library", .libPaths()))
library(readxl)

file_path <- "data/raw/Advanced_lines-2021 Orignal .xlsx"
column_names <- c(
  "entry", "plant_height", "tillers", "panicle_length", "filled_grains",
  "unfilled_grains", "grain_weight_100", "plant_weight", "recorded_yield"
)
measurement_columns <- column_names[-1]
replications <- c("R-1", "R-2", "R-3")

# A function lets us repeat the same import steps for each worksheet.
import_replication <- function(sheet_name) {
  # trim_ws = FALSE keeps spaces in the original entry text.
  # The first column is text; the other eight columns are measurements.
  imported <- read_excel(
    path = file_path, sheet = sheet_name, range = "A3:I55",
    col_names = FALSE, col_types = c("text", rep("numeric", 8)),
    trim_ws = FALSE
  )
  # A base R data.frame prints every flagged record without tibble truncation.
  imported <- as.data.frame(imported)
  names(imported) <- column_names
  imported$entry_original <- imported$entry
  # trimws() removes surrounding whitespace, leaving internal text unchanged.
  imported$entry <- trimws(imported$entry_original)
  imported$replication <- sheet_name
  imported
}

rice_r1 <- import_replication("R-1")
rice_r2 <- import_replication("R-2")
rice_r3 <- import_replication("R-3")

# rbind() stacks tables with matching columns underneath one another.
phenotype_2021 <- rbind(rice_r1, rice_r2, rice_r3)
cat("\nTotal rows:", nrow(phenotype_2021), "\n")
cat("Unique nonmissing entries:",
    length(unique(phenotype_2021$entry[!is.na(phenotype_2021$entry)])), "\n")
cat("Missing entry labels:", sum(is.na(phenotype_2021$entry)), "\n")

# table() counts how many records have each replication label.
cat("\nRecords per replication:\n")
print(table(phenotype_2021$replication, useNA = "ifany"))

# paste0() builds the expected labels E-1, E-2, ..., E-53.
# setdiff(x, y) shows values in x that are absent from y.
expected_entries <- paste0("E-", 1:53)
cat("\nExpected entry check for each replication:\n")
for (replication_name in replications) {
  entries <- phenotype_2021$entry[
    phenotype_2021$replication == replication_name
  ]
  missing_entries <- setdiff(expected_entries, entries)
  unexpected_entries <- setdiff(entries, expected_entries)
  cat("\n", replication_name, "\n", sep = "")
  cat("Contains all E-1 through E-53:", length(missing_entries) == 0, "\n")
  cat("Missing expected entries:\n")
  print(missing_entries)
  cat("Unexpected entry labels (including NA, if present):\n")
  print(unexpected_entries)
}

# duplicated() flags repeated combinations after their first occurrence.
entry_keys <- phenotype_2021[c("entry", "replication")]
duplicate_later <- duplicated(entry_keys)
# Checking from both ends includes the first record of each duplicate group.
duplicate_any <- duplicate_later | duplicated(entry_keys, fromLast = TRUE)
cat("\nDuplicate entry-by-replication combinations beyond the first record:",
    sum(duplicate_later), "\n")
cat("All records belonging to duplicate combinations:\n")
print(phenotype_2021[duplicate_any, ], row.names = FALSE)

# is.na() gives TRUE for missing values. sum() counts TRUE values as 1.
# na.rm = TRUE excludes missing values when counting zeros.
measurements <- phenotype_2021[measurement_columns]
measurement_audit <- data.frame(
  measurement = measurement_columns,
  missing = sapply(measurements, function(x) sum(is.na(x))),
  zeros = sapply(measurements, function(x) sum(x == 0, na.rm = TRUE)),
  row.names = NULL
)
cat("\nMissing values and zeros by measurement:\n")
print(measurement_audit, row.names = FALSE)

# rowSums() counts flags across the eight measurements for each record.
has_missing <- rowSums(is.na(measurements)) > 0
has_zero <- rowSums(measurements == 0, na.rm = TRUE) > 0
cat("\nRecords containing missing measurements or zeros:",
    sum(has_missing | has_zero), "\n")
print(phenotype_2021[has_missing | has_zero, ], row.names = FALSE)

# Save every imported row, including flagged observations. No values are fixed.
# entry_original retains the entry text before trimws(). NA is written as NA.
dir.create("outputs", showWarnings = FALSE)
write.csv(phenotype_2021, "outputs/phenotype_2021_imported.csv",
          row.names = FALSE, na = "NA")
cat("\nSaved outputs/phenotype_2021_imported.csv\n")
