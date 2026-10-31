# Reproducibility verification

Verification used Rscript --vanilla run_all.R in a fresh R process from the project root. No package installation occurred. Every analysis script 00 through 12 completed successfully in a separate source environment; no prior-session objects were supplied.

Pre-run saved CSV tables compared: 61
Equal within numerical tolerance 1e-10: 61
Exactly identical parsed tables: 61
Remaining numerical comparison failures: 0
Raw input MD5 checksums unchanged versus pre-run snapshot: TRUE

The pre-run snapshot was taken before regeneration, covering every existing output CSV. Comparison checks columns, row order, values, and attributes after parsing CSV; numerical equality uses all.equal tolerance 1e-10. Exact parsed equality is reported separately. New session/checksum/comparison files have no pre-run counterpart. No numerical differences were found requiring investigation.
This verifies reproducibility with the current R/readxl environment and input files. It does not establish scientific validity, portability across every R/package version, or byte-identical PNG/RDS serialization. Prior visual inspection is documented in the step notes; this verification compares numerical CSVs, not pixel hashes.

MD5 checksums are change-detection records, not authentication guarantees. outputs/input_file_checksums.csv includes both raw workbooks and their preserved root copies; outputs/analysis_script_checksums.csv records runner and script versions. A separate comparison result is in outputs/reproducibility_comparison.csv. outputs/sessionInfo.txt records the pipeline session.

Remaining execution failures: none. Scientific unresolved issues remain as described in the README and results brief (units, sampling/design provenance, yield discrepancies, residual variance, candidate selection optimism, and E-42/E-45 provenance).

## Pipeline session

```
R version 4.5.2 (2025-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 26.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.1 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.1;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
 [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
 [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
[10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   

time zone: Asia/Karachi
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] readxl_1.6.0

loaded via a namespace (and not attached):
 [1] compiler_4.5.2   magrittr_2.0.5   cli_3.6.6        tools_4.5.2     
 [5] pillar_1.11.1    glue_1.8.1       tibble_3.3.1     utf8_1.2.6      
 [9] cellranger_1.1.0 vctrs_0.7.3      lifecycle_1.0.5  rematch_2.0.0   
[13] pkgconfig_2.0.3  rlang_1.3.0     
```
