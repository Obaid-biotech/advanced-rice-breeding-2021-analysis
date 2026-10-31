# 2021 yield consistency audit

Run from the project folder with `Rscript R/02_yield_consistency.R`.
The script uses the existing `.R-library`, readxl, and base R. It checks the
imported CSV's 11 column names and reads Combined cells A3:E55 from the original
workbook. Source data and earlier scripts are preserved.

The factor 135.85 is an **unverified numerical relationship**, not a validated
conversion to yield per hectare. Differences are flagged only when their
absolute value exceeds 0.000001. The factor difference is recorded yield minus
135.85 times plant weight.

| Entry | Replication | Plant weight | Recorded yield | Expected from factor | Difference (rounded) |
|---|---|---:|---:|---:|---:|
| E-31 | R-1 | 28.3466666666667 | 4850.89466666667 | 3850.89466666667 | 1000 |
| E-39 | R-3 | 34.3866666666667 | 5671.42866666667 | 4671.42866666667 | 1000 |
| E-40 | R-3 | 62.7 | 5517.795 | 8517.795 | -3000 |

All 159 factor differences were checkable; three exceeded tolerance.
E-46 in R-3 has plant weight 0 and recorded yield 0. This satisfies the equation
because 135.85 times 0 is 0. Arithmetic agreement cannot establish what the zeros
mean biologically. Field/source records are needed to clarify whether these are
observations, missing-data codes, or recording issues. The zeros remain intact.

After trimming whitespace, Combined has 53 unique, nonmissing IDs, exactly E-1
through E-53. There are no blank or duplicate IDs. Matching by entry and
replication gives 159 matched keys and no unmatched keys. Two yields differ:

| Entry | Replication | Combined yield | Replication yield | Combined minus replication |
|---|---|---:|---:|---:|
| E-42 | R-2 | 6322.171 | 8322.171 | -2000 |
| E-46 | R-3 | 130 | 0 | 130 |

All 53 Combined means equal the arithmetic mean of their own three yield columns
within tolerance; none is uncheckable. This checks internal arithmetic and does
not resolve the differences between sheets.

Reports in `outputs/`:

- `yield_factor_discrepancies.csv`: three factor flags.
- `yield_zero_records.csv`: the one record with a zero measurement in this check.
- `yield_between_sheet_discrepancies.csv`: two yield differences; also includes
  unmatched keys or missing yields if encountered on a future run.
- `yield_combined_mean_discrepancies.csv`: no flagged rows (headers only).

`abs()` returns magnitude: `abs(-3000)` is 3000. The comparison `> 0.000001`
flags a difference outside tolerance; `== 0` tests for zero. `&` requires both
conditions, while `|` requires at least one. For example,
`phenotype[factor_flag, report_columns]` selects rows where `factor_flag` is TRUE
and displays the six requested columns. Missing differences are counted
separately, so they are not silently treated as passing the factor check.

At this historical yield-audit stage (step 02), no values were corrected,
zeros replaced, lines ranked, or models fitted. Provisional phenotype modeling
was subsequently performed in steps 07–09, followed by exploratory candidate
selection in step 10. The audit findings above remain unchanged.
