# Candidate marker profiles

Run `Rscript R/11_candidate_marker_profiles.R` (base R). Existing CSV structures were inspected first. All four top-ten membership sets from script 10 agree; all ten candidates, including E-37, remain included. Mapping is validated by explicit IDs, confirmed suffixes, and sample names. Correspondence remains user-confirmed, not independently field-verified.

## Definitions

Panels M10, M20, and M30 retain markers with missing calls / 53 original samples <= 0.10, 0.20, or 0.30, respectively, and exactly two observed alleles across the same original 53 samples. These are exploratory panels. Missingness and allele counts are recalculated from the audited calls and checked against audit 04. Sample filtering is never applied. The membership CSV includes every original marker for every panel.
Original genotype calls are read alongside their documented canonical representation and revalidated. Valid A/C/G/T pairs are alphabetically sorted: G:A and A:G are equivalent. A:A and A:G remain distinct genotypes. ?, Uncallable, Bad, and blanks remain missing. Unexpected strings stop the script. No calls are imputed.
For each of 45 unordered candidate pairs and each panel, the denominator is the number of panel markers with valid calls in BOTH samples. Mismatch proportion = differing canonical genotypes / jointly called markers. No joint calls yields NA. The available fraction is jointly called markers / all panel markers, not / 502. Different pairs and panels can use different marker subsets.
Limited-data review flags identify denominators or available fractions at or below the panel-specific lower quartile among the 45 pairs, plus any pair with no joint calls. The quartile reference values and continuous counts/fractions are exported. This is an explicitly relative descriptive flag, NOT a validated reliability cutoff; unflagged comparisons are not certified reliable.
Candidate QC metrics and flags for exceeding the earlier exploratory 10%, 20%, and 30% sample-missingness limits are preserved in a separate CSV. No candidate is removed because of QC. Missingness uses 502 original markers; heterozygosity uses valid calls as denominator.

## Actual results

 panel panel_markers pairs min_joint_calls max_joint_calls
   M10           364    45             246             364
   M20           432    45             266             430
   M30           446    45             271             443
 min_fraction_available max_fraction_available min_mismatch max_mismatch
              0.6758242              1.0000000            0    0.6074074
              0.6157407              0.9953704            0    0.6250000
              0.6076233              0.9932735            0    0.6187291
 relative_review_flags no_joint_calls
                    12              0
                    12              0
                    12              0

Largest across-panel changes (not a ranking of crosses):

 entry_1 entry_2 M10_mismatch_proportion M20_mismatch_proportion
     E-3    E-45               0.2670455               0.3000000
     E-3    E-42               0.2691218               0.3017032
     E-3    E-10               0.2096317               0.2432432
     E-3    E-11               0.5738636               0.6063570
     E-3    E-51               0.2215909               0.2536585
 M30_mismatch_proportion panel_mismatch_range
               0.3080569           0.04101142
               0.3096927           0.04057086
               0.2457831           0.03615140
               0.6076555           0.03379187
               0.2517815           0.03206763

## Interpretation and outputs

All three heatmaps were visually inspected. Broad patterns persist across panels, with E-11 showing higher mismatch against other candidates than most other displayed pairs. E-42/E-45 has zero observed mismatch across 363, 430, and 443 jointly called markers in M10, M20, and M30. This does not establish genome-wide identity, duplicate sample status, or pedigree. The largest across-panel range is about 0.0410 for E-3/E-45 (4.10 percentage points).
Each panel has 12 relative lower-quartile review flags, including all nine pairs involving E-37. E-37 remains retained, with its original QC metrics and exploratory flags exported. A lower available fraction means a larger part of that panel could not contribute to the pair comparison; a zero mismatch is agreement among observed calls only.
These quantities describe observed genotype mismatch in a targeted marker panel. They are not genome-wide genetic distances, causal marker effects, favorable alleles, or evidence that a cross will produce heterosis. Panel selection and missing data affect which loci contribute, so variation across panels is descriptive sensitivity evidence. No crosses are recommended from mismatch alone.
All three heatmaps use the same fixed numeric entry order and the same 0–1 color scale, without clustering. Cell labels show mismatch proportions rounded to two decimals. Diagonals are not pair comparisons. Joint-call matrices and long reports provide exact denominators and available fractions.
Outputs under outputs/: candidate_marker_panel_membership.csv; candidate_marker_pairwise_profiles.csv (135 rows); candidate_marker_panel_comparison.csv (45 rows); candidate_marker_panel_summary.csv; candidate_marker_profiles_sample_qc.csv; and M10/M20/M30 mismatch matrices, joint-count matrices, and heatmap PNGs with the candidate_marker_ prefix.
Source observations and prior outputs are preserved. No candidates are removed, calls imputed, associations tested, or crosses recommended.

R version 4.5.2 (2025-10-31)
