# 2021 exploratory results brief

## Dataset and quality audit

The 2021 phenotype workbook contributes 53 entries in three reported RCBD replications (159 records). Imported measurements have no missing values or duplicate entry-by-replication keys. Original ID text is retained. E-46/R-3's zero grain weight per plant remains unchanged, without declaring it crop failure or missing data.

The user confirms `plant_weight` means grain weight per plant and E-n corresponds to NYUT_21_n for all 53 entries. Units and sampling details remain unresolved; sample correspondence has not been independently field-verified.

`recorded_yield` remains separate. The relationship 135.85 × plant_weight is unverified, not a hectare conversion. Three records differ from it: E-31/R-1 and E-39/R-3 by approximately +1000, E-40/R-3 by −3000. Combined differs from replication records for E-42/R-2 by −2000 and E-46/R-3 by +130. Its own three-column means are internally consistent. These checks locate discrepancies without correcting them.

## Provisional phenotype models

Both response scales use additive fixed effects: response ~ replication + entry. Scenario A includes all observations; B omits only E-46/R-3 as explicitly requested sensitivity. All four models have rank 55; residual df are 104 in A and 103 in B.

| Scale | Scenario | Block F; p | Entry F; p |
|---|---|---|---|
| Original | A | 8.238; 0.000477 | 4.389; 7.27e-11 |
| Original | B | 7.997; 0.000591 | 4.024; 8.78e-10 |
| Square root | A | 9.912; 0.000115 | 5.885; 9.62e-15 |
| Square root | B | 8.822; 0.000291 | 4.812; 5.85e-12 |

Partial tests use drop1(), retaining the other term, rather than order-dependent sequential tests. Small p-values are conditional on model assumptions and do not establish genetic superiority. Residual-versus-fitted and scale-location diagnostics show increasing spread on the original scale. Square-root diagnostics reduce that pattern, with remaining departures; E-46/R-3 becomes more unusual on the transformed scale. No assumption is declared solved. Diagnostics cannot verify randomization or independence.

Entry estimates average equally across all three blocks, including model prediction for E-46's omitted block in B. Full coefficient covariances are used for standard errors. Original A means agree with raw entry means within 1.4e-13. Only E-46 changes beyond rounding between scenarios: approximately +0.1252 original units and +0.6201 square-root units. Its uncertainty increases in B. Pointwise t intervals are for means, not prediction intervals, simultaneous coverage, selection-adjusted inference, or pairwise tests. Original E-46 intervals extend below zero because the unconstrained model intervals are symmetric; they are not physical claims.

Square-root means target a different quantity and remain on their modeled scale. Squaring a fitted transformed mean omits a variance component and does not directly recover arithmetic mean grain weight. Residual SD and ordinary AIC must not be compared directly across response scales. No model is selected by smaller p-values.

## Candidate membership sensitivity

Completed selection is exploratory **grain-weight-only** prioritization. Sizes 5 and 10 are explicit choices, not biologically validated cutoffs. An absolute 1e-10 cutoff-tie tolerance yields actual sizes 5 and 10 in all models.

- Five: E-3, E-15, E-40, E-45, E-51.
- Ten: E-3, E-10, E-11, E-15, E-16, E-37, E-40, E-42, E-45, E-51.

Membership agrees across both scales and scenarios (5/5 and 10/10 overlaps). This is sensitivity evidence from closely related fits, not independent validation or a probability of superiority. Candidate selection on the same data may exaggerate performance. Other recorded traits have descriptive means only: **multi-trait selection was not performed**, and no favorable directions, weights or composite objective were assigned. The evidence table includes specific yield and molecular QC cautions for each candidate.

## KASP characterization

The original call region has 53 samples × 502 markers (26,606 cells). Missing labels account for 2,039 calls (7.66%): 2,003 ?, 32 Uncallable and four Bad. No blank or unexpected call strings were found. Canonical allele pairs equate reversed order while preserving homozygous/heterozygous distinctions. Original calls remain saved.

There are 48 monomorphic and 454 biallelic markers, no multiallelic or entirely missing markers, and no missing/duplicate marker names, missing positions, or repeated positions. Syntax checks do not verify a reference genome build.

Sample missingness ranges 1.59–30.88%; heterozygosity among valid calls ranges 0.22–17.78%. E-37 has 30.88% missingness and remains retained. Earlier nine sample/marker QC scenarios retain 42–52 samples and 381–450 polymorphic biallelic markers, with no final scenario chosen.

For the ten candidates, exploratory marker panels use missingness across all 53 original samples plus two observed alleles. At 10/20/30% limits, 364/432/446 markers qualify. Pairwise jointly called counts range 246–364, 266–430 and 271–443, respectively. Mismatch ranges are 0–60.74%, 0–62.50%, and 0–61.87%. Each panel flags 12 pairs in its relative lower-quartile shared-data review group, including all nine E-37 pairs; this is not a validated reliability cutoff. Largest panel mismatch change: E-3/E-45, about 4.10 percentage points.

These are observed genotype mismatches in targeted panels, not genome-wide genetic distances, causal effects, favorable alleles or evidence for heterosis. No candidates are removed or calls imputed.

## E-42/E-45 investigation

All 502 original markers are classified: 487 matching valid calls, zero mismatching valid calls, two missing only in E-42, one only in E-45, and 12 missing in both. The disjoint categories sum to 502. Laboratory labels are RC21075 (B06) and RC21080 (E06), both plate TSK14060. Original phenotype records differ and are exported side by side, keeping grain weight distinct from recorded_yield.

Neither label/well differences nor observed marker agreement independently establishes biological identity. Possible explanations include distinct materials insufficiently distinguished by this panel, shared ancestry, handling/labelling issues, or assay/export issues. Review field and seed records, pedigrees/accessions, collection chain of custody, extraction logs, plate maps, raw fluorescence/cluster files and export history. Independently sampled traceable material, repeated extractions with controls, additional informative markers or sequencing, and within-entry replicate sampling can help distinguish alternatives. No identity is corrected and no crosses recommended.

## Conclusions and unresolved questions

The completed workflow supports review of 2021 grain-weight candidate sensitivity and descriptive molecular evidence, with explicit uncertainty. It does not establish multi-year stability, broad adaptation, validated yield per hectare, biological identity equivalence, genetic superiority, or causal marker–trait links.

Next steps are to resolve units, sampling and discrepancy provenance; confirm the experimental design; investigate laboratory/seed identity records; define a breeding objective for any multi-trait work; and obtain independent trials. Model assumptions, selection optimism and targeted-panel ascertainment constrain interpretation. No marker association tests, significance letters, or cross recommendations have been made.

Reproducibility verification and remaining execution failures are documented in `notes/13_reproducibility.md`. Session information and input checksums are saved under `outputs/`.
