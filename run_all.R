# Run from the project root: Rscript --vanilla run_all.R
# Install readxl separately (see README). Analysis never installs packages.
local({
  if (!file.exists("R/00_first_import.R") || !dir.exists("data/raw")) {
    stop("Run from the Advanced_Lines_Project root directory.")
  }
  .libPaths(c(".R-library", .libPaths()))
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl separately; see README.md.")
  scripts <- sprintf("R/%s", c(
    "00_first_import.R", "01_phenotype_audit.R", "02_yield_consistency.R",
    "03_sample_mapping.R", "04_kasp_audit.R", "05_kasp_qc_sensitivity.R",
    "06_phenotype_exploration.R", "07_grain_weight_rcbd.R",
    "08_sqrt_model_sensitivity.R", "09_entry_estimates.R",
    "10_candidate_sensitivity.R", "11_candidate_marker_profiles.R",
    "12_candidate_evidence_summary.R"))
  stopifnot(all(file.exists(scripts)))
  inputs <- c(list.files("data/raw", pattern = "\\.xlsx$", full.names = TRUE),
              list.files(".", pattern = "\\.xlsx$", full.names = TRUE))
  before <- tools::md5sum(inputs)
  for (script in scripts) {
    cat("\nRunning", script, "\n")
    # Isolate script variables so one script cannot overwrite another's state.
    source(script, local = new.env(parent = globalenv()), echo = FALSE)
  }
  after <- tools::md5sum(inputs)
  stopifnot(identical(before, after))
  write.csv(data.frame(file = names(after), md5 = unname(after)),
            "outputs/input_file_checksums.csv", row.names = FALSE)
  write.csv(data.frame(file = c("run_all.R", scripts),
                      md5 = unname(tools::md5sum(c("run_all.R", scripts)))),
            "outputs/analysis_script_checksums.csv", row.names = FALSE)
  writeLines(capture.output(sessionInfo()), "outputs/sessionInfo.txt")
  cat("\nAll scripts completed; input workbook checksums unchanged.\n")
})
