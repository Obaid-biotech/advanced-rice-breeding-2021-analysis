# Data dictionary

## Sample identifiers

`outputs/sample_id_mapping.csv` contains one row per entry (53 rows).

| Column | Meaning | Source |
|---|---|---|
| phenotype_entry | Phenotype entry ID, E-1 through E-53 | Distinct `entry` values in `outputs/phenotype_2021_imported.csv` |
| kasp_id | KASP NIBG identifier, NYUT_21_1 through NYUT_21_53 | `data/raw/Advanced_lines-2021-KASP_GenotypeData.xlsx`, Sheet1, A3:A55 |
| sample_name | Exact text in the KASP Sample Name column, such as RF21001 or RC21099 | Same workbook, Sheet1, B3:B55 |

On 2026-10-30, the user confirmed that phenotype E-n corresponds to KASP
NYUT_21_n for all 53 entries. This correspondence is **user-confirmed**;
it has not been independently verified from field records. Preserve the exact
KASP prefix `Advanced_lines_21_`, even though other project names use NUYT.

`entry_original` in the phenotype import preserves the original entry text;
`entry` is its whitespace-trimmed working ID. The three replication records
per entry refer to the same mapped sample; they are not three separate KASP IDs.
Sample names are taken from the workbook, not inferred from entry numbers.

## Phenotype meaning and interpretation

The user confirms `plant_weight` means grain weight per plant. Its unit and
sampling details remain unconfirmed. `recorded_yield` remains a separate recorded
measurement; the numerical relationship 135.85 times `plant_weight` is unverified
and is not a validated conversion to yield per hectare. Audit 02 reports factor
and between-sheet discrepancies without correcting them. The E-46/R-3 zero is
preserved without declaring biological failure or missingness.

Original-response entry estimates concern mean grain weight equally averaged
over three trial blocks. Square-root estimates concern mean square-root grain
weight and are retained on that modeled scale. Their pointwise confidence
intervals are provisional and model-dependent.
