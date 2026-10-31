# Phenotype descriptive exploration

Run `Rscript R/06_phenotype_exploration.R` from the project folder. Only base R is used. The expected 11 columns and eight numeric measurement columns are validated before summarizing.

## Meaning and unresolved questions

The user confirms that plant_weight means grain weight per plant. Its physical unit remains unconfirmed, so both plots label it 'unit unconfirmed'. Original column names and data are preserved.
E-46/R-3 has a recorded zero. It is retained in summaries and plots without interpreting it as confirmed crop failure or as missing data. Its biological meaning requires clarification.
recorded_yield is summarized separately, never substituted for grain weight per plant. The factor 135.85 is an unverified numerical relationship, not a validated conversion to yield per hectare. Audit 02 found factor discrepancies for E-31/R-1 (+1000), E-39/R-3 (+1000), and E-40/R-3 (-3000), rounded. Combined differs from replication recorded_yield for E-42/R-2 (-2000) and E-46/R-3 (+130), expressed as Combined minus replication. See notes/02_yield_consistency.md. No correction is made here.

## Grain weight per plant summaries

| Replication | Observed | Missing | Zeros | Mean | SD | Minimum | Median | Maximum |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| R-1 | 53 | 0 | 0 | 35.412 | 11.653 | 3.830 | 34.567 | 74.340 |
| R-2 | 53 | 0 | 0 | 37.160 | 10.842 | 6.920 | 36.120 | 61.260 |
| R-3 | 53 | 0 | 1 | 31.280 | 11.025 | 0.000 | 29.790 | 62.700 |

## How to read the plots

The replication plot shows every individual observation, with horizontal jitter to reduce overlap. The vertical values are unchanged. Boxes show the median and middle half of the data; whiskers extend to the most extreme observations within 1.5 times the interquartile range. All observations, including those beyond whiskers, appear as points. The E-46/R-3 zero is highlighted.
The profile plot highlights E-46: R-1 = 3.830; R-2 = 6.920; R-3 = 0.000
Other entries are subdued but all are shown. Connected profiles reveal within-entry variation across the three replications. Replication is a grouping label, not a time axis. Crossing lines show that entries do not follow identical profiles.
Differences between the displayed means, spreads, and profiles are descriptive observations only. These plots establish neither statistically significant replication differences nor genetic superiority. No lines are ranked.

## Functions and denominators

Grouping uses phenotype[phenotype$replication == replication_name, ] to select one replication's rows. Each measurement is then summarized separately within that group.
mean() adds the observed values and divides by their count. sd() measures spread around the mean, using the sample standard deviation with denominator n - 1. median() returns the middle sorted value (or the average of the two middle values for an even count).
Missing values are counted separately and omitted from these calculations; zeros are observed values and remain included. SD is NA with fewer than two observations; other summaries are NA when there are no observations. The output records number_observed, missing_count, and zero_count explicitly.

## Outputs

- outputs/phenotype_summary_by_replication.csv: all eight measurements, each summarized in all three replications.
- outputs/grain_weight_per_plant_summary_by_replication.csv: the plant_weight subset.
- outputs/grain_weight_per_plant_by_replication.png: distributions and individual observations.
- outputs/grain_weight_per_plant_profiles.png: one connected profile per entry.
No observations were removed or replaced, outliers filtered, or statistical models fitted. All previous files remain unchanged.
