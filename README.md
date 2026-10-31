# Advanced 2021 exploratory grain-weight and KASP analysis

## Research question and completed objectives

Within the reported 2021 trial, which entries merit exploratory follow-up based on grain weight per plant, how sensitive is membership to response scale and the E-46/R-3 zero, and what do the available targeted KASP calls describe about those candidates?

Completed work covers phenotype and yield-consistency audits, a sample-ID crosswalk, KASP characterization and QC sensitivity, descriptive phenotype exploration, provisional fixed-effects RCBD models, equal-block entry means and pointwise intervals, exploratory five- and ten-entry selection, targeted-panel mismatch profiles, and a candidate evidence summary. Grain-weight selection is complete as an exploratory exercise. **Multi-trait selection has not been performed**: other traits are described without favorable directions, weights, or a composite score.

## Inputs and correspondence

- `data/raw/Advanced_Lines-2021 Orignal .xlsx`: R-1/R-2/R-3, A3:I55; Combined, A3:E55.
- `data/raw/Advanced_Lines-2021-KASP_GenotypeData.xlsx`: Sheet1 sample rows 3–55; marker names in row 1, positions in row 2, marker columns E:SL (5–506). A:D contain laboratory metadata. Summary cells are excluded.

The original workbooks, including root-level copies, are preserved. The user confirms phenotype E-n corresponds to **Advanced_lines_21_n** for all 53 entries. This is user-confirmed correspondence, not independent verification from field records. The exact KASP spelling is retained. See `notes/data_dictionary.md` and `notes/decision_notes.md`.

The user confirms `plant_weight` means grain weight per plant.

## Software and setup

Verified environment: R 4.5.2 and readxl 1.6.0. Only readxl and base/recommended R facilities are used. The existing `.R-library` is loaded by import scripts and the runner. Package installation is separate from analysis; if readxl is absent, run once from the project root:

```r
dir.create(".R-library", showWarnings = FALSE)
install.packages("readxl", lib = ".R-library", repos = "https://cloud.r-project.org")
```

Then, from `/mnt/d/Advanced_lines_Project`:

```sh
Rscript --vanilla run_all.R
```

This regenerates analysis outputs and generated step notes, stops on errors, checks raw workbooks remain unchanged, and writes `outputs/sessionInfo.txt`, `input_file_checksums.csv` (MD5), and `analysis_script_checksums.csv`. It never installs packages. Script 12 uses a deliberately narrow base-R ZIP/XML reader for the inspected workbook's A:D shared-string metadata; a different Excel layout may require revision. No current-session objects are required.

## Execution order and principal outputs

| Script | Purpose and main outputs in `outputs/` |
|---|---|
| 00 | First-import lesson; prints workbook inspection |
| 01 | `phenotype_2021_imported.csv`; uncorrected import with original IDs |
| 02 | `yield_*discrepancies.csv`, `yield_zero_records.csv` |
| 03 | `sample_id_mapping.csv` |
| 04 | `kasp_*`: original/canonical calls, sample/marker QC, metadata and allele counts |
| 05 | `kasp_qc_sensitivity.csv`, marker details and PNG |
| 06 | Phenotype summaries and grain-weight distribution/profile PNGs |
| 07 | `grain_weight_rcbd_*`: original-response models' partial tests and diagnostics |
| 08 | `sqrt_sensitivity_*`: both response scales, diagnostics, tests and model RDS copies |
| 09 | `entry_estimates_*`: equal-block means, CIs, design rows and validation |
| 10 | `candidate_sensitivity_*`, membership, ordered estimates and overlaps |
| 11 | `candidate_marker_*`: exploratory panel membership, pair mismatch and denominators, heatmaps |
| 12 | `candidate_evidence_summary.csv`, `e42_e45_*` marker, metadata and phenotype tables |

Each script's full filename is in `R/`; execute 00 through 12 in numeric order. Dependency inspection: 01 produces the phenotype input used by 02/03/06–09/12; 03 feeds 04/10–12; 04 feeds 05/10–12; 09 feeds 10; 10 feeds 11/12; 11 supplies candidate QC to 12; 02 supplies yield flags to 12. Scripts 07/08 inform 09's definitions but 09 refits rather than requiring session models. Separate environments in the runner prevent script-variable collisions. Step notes are under `notes/`.

## Findings from saved results

- Phenotypes: 159 records, 53 entries, three complete replications; no missing measurements or duplicate entry/block keys. E-46/R-3 has recorded zero grain weight and recorded yield; neither is corrected or given a confirmed biological interpretation.
- Yield audit: factor discrepancies for E-31/R-1 (+1000), E-39/R-3 (+1000), E-40/R-3 (−3000), rounded. Combined minus replication differences: E-42/R-2 −2000 and E-46/R-3 +130. Combined means pass their internal arithmetic check.
- KASP: 53 × 502 calls, 2,039 missing (7.66%), no unexpected strings; 48 monomorphic and 454 biallelic markers; no missing/duplicate marker names or repeated positions.
- Models: `response ~ replication + entry`, both factors. A uses 159 observations; B omits only E-46/R-3 as sensitivity. Rank is 55; residual df are 104/103. Original residual SDs are 7.660/7.697. Residual spread increases with fitted values. Square roots reduce that pattern but do not establish valid assumptions.
- Five-entry membership in all four models: E-3, E-15, E-40, E-45, E-51. Ten-entry membership adds E-10, E-11, E-16, E-37, E-42. No cutoff ties expand the sets at absolute tolerance 1e-10. This is exploratory prioritization, not statistical superiority.
- Candidate panels at 10/20/30% missingness over all 53 samples contain 364/432/446 biallelic markers. Shared-call denominators are exported for every pair; E-37 remains included with its 30.88% sample missingness flag.
- E-42/E-45 match at all 487 jointly valid original markers; no mismatches, two missing only in E-42, one only in E-45, 12 in both. Their labels are RC21075/B06 and RC21080/E06 on plate TSK14060. These facts do not establish identical biological identity.

See `notes/13_results_brief.md` for interpretation and `notes/13_reproducibility.md` for the fresh-process comparison against saved results.

## Limitations and follow-up

Only 2021 is analyzed; stable performance across years or broad adaptation is not established. Grain-weight units, experimental-unit/sampling details, randomization and independence need field-record confirmation. Residual variance concerns make conventional F-tests and pointwise intervals provisional. Overlapping intervals are not pairwise tests. Square-root means remain on that scale; squaring them does not recover arithmetic means. Residual SD and ordinary AIC are not directly comparable across response scales, and model choice must not chase smaller p-values.

Candidate agreement is sensitivity evidence from related fits, not independent validation or a probability of superiority. Selection can exaggerate performance. No multi-trait objective is defined. KASP mismatch describes observed calls in targeted panels, not genome-wide distance, favorable alleles, causal effects or heterosis. Relative lower-quartile shared-data flags are review aids, not validated reliability cutoffs.

Follow-up: resolve units, sampling and yield discrepancies; review field, seed, accession, collection and laboratory provenance; repeat traceable DNA sampling/genotyping and use additional informative markers to investigate E-42/E-45 alternatives; confirm performance in independent trials and environments; define a breeding objective before multi-trait selection. No identities are corrected and no crosses recommended from mismatch alone. No associations, significance letters or heritability estimates were produced.
