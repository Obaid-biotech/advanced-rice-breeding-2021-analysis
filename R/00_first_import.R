# First lesson: import only the 2021 rice data, without cleaning it.
# Run from the project folder: Rscript R/00_first_import.R

# Use the project's user-owned package library, if it exists.
# .libPaths() tells R where to look for installed packages.
if (dir.exists(".R-library")) {
  .libPaths(c(".R-library", .libPaths()))
}

# library() loads the package that provides Excel-reading functions.
# If needed, run this installation command once from the project folder:
# install.packages("readxl", lib = ".R-library",
#                  repos = "https://cloud.r-project.org")
library(readxl)

# <- assigns the value on its right to the name on its left.
# This relative path starts from the project's working directory.
file_path <- "data/raw/Advanced_lines-2021 Orignal .xlsx"

# Show all worksheet names before importing one worksheet.
cat("Worksheet names:\n")
print(excel_sheets(file_path))

# Read only R-1, cells A3:I55. FALSE means the first row is data,
# not column headings. Reading the workbook does not modify it.
rice_2021 <- read_excel(
  path = file_path,
  sheet = "R-1",
  range = "A3:I55",
  col_names = FALSE
)

# Give the nine imported columns the agreed names.
names(rice_2021) <- c(
  "entry", "plant_height", "tillers", "panicle_length", "filled_grains",
  "unfilled_grains", "grain_weight_100", "plant_weight", "recorded_yield"
)

# head() displays the first six rows by default.
cat("\nFirst six rows:\n")
print(head(rice_2021))

# dim() gives the number of rows followed by the number of columns.
cat("\nDimensions (rows, columns):\n")
print(dim(rice_2021))

cat("\nColumn names:\n")
print(names(rice_2021))

# str() shows the structure, including each column's type and sample values.
cat("\nStructure and column types:\n")
str(rice_2021)

# $ selects a named column. This example is a comment, so it is not run:
# head(rice_2021$plant_height)
