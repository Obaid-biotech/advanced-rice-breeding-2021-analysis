# 2021 KASP audit

Run `Rscript R/04_kasp_audit.R` from the project folder.

Sheet1 rows 3–55 and columns 5–506 (E:SL) contain 53 samples and 502 markers.
Marker names come from row 1 and positions from row 2. Headers, sample metadata, and summary cells are excluded from the call matrix.
Sample IDs and Sample Name text agree with sample_id_mapping.csv: TRUE
The identity correspondence remains user-confirmed, not independently verified from field records.

## Results

samples: 53
markers: 502
total_cells: 26606
missing_calls: 2039
unexpected_calls: 0
monomorphic_markers: 48
biallelic_markers: 454
multiallelic_markers: 0
entirely_missing_markers: 0
markers_without_valid_calls: 0
missing_marker_names: 0
markers_with_duplicate_names: 0
missing_positions: 0
unexpected_position_formats: 0
markers_with_repeated_positions: 0
distinct_repeated_positions: 0

Sample missing proportion range: 0.0159362549800797 to 0.308764940239044
Marker missing proportion range: 0 to 0.981132075471698
Sample heterozygous proportion range: 0.00224719101123596 to 0.177824267782427

Every original call string and its frequency:

?: 2003
A:A: 5886
A:C: 3
A:G: 2
A:T: 1
Bad: 4
C:A: 78
C:C: 6039
C:T: 6
G:A: 267
G:C: 68
G:G: 6366
G:T: 2
T:A: 71
T:C: 189
T:G: 48
T:T: 5541
Uncallable: 32

## Definitions and denominators

Missing calls are ?, Uncallable, Bad, and blank cells (including whitespace-only cells). Sample missing proportion divides missing cells by all 502 markers. Marker missing proportion divides missing cells by all 53 samples.
Valid calls are exact A/C/G/T pairs separated by a colon. Unexpected strings are reported separately, not silently counted as the defined missing labels; they are unscorable in the derived copy.
Sample heterozygous proportion divides valid pairs with two different alleles by all valid pairs for that sample. Missing and unexpected cells are excluded from this denominator. If there are no valid calls, the proportion is NA.
Observed allele counts use distinct A/C/G/T letters across valid pairs only. One allele means monomorphic, two biallelic, and more than two multiallelic. Entirely missing means all 53 cells have a defined missing label; markers with no valid calls are also counted separately.
The derived analysis copy sorts valid pairs alphabetically, so G:A and A:G both become A:G. Original text is retained alongside it. Repeated positions are reported without removing markers. Position syntax checks do not verify a reference genome or genome build.

## Outputs

All audit reports have the outputs/kasp_ prefix: sample_id_validation, marker_metadata_audit, repeated_positions, original_call_frequencies, unexpected_calls, sample_audit, marker_audit, original_and_analysis_calls, audit_summary, and software_versions (all CSV). Empty reports retain column headers.
The original_and_analysis_calls report retains every cell and its status. Original blank cells are written as NA; their original cell locations are retained through sample IDs and Excel marker columns.
No samples or markers were filtered, calls imputed, marker calls joined to phenotype measurements, or marker–trait associations tested.

R: R version 4.5.2 (2025-10-31)
readxl: 1.6.0
