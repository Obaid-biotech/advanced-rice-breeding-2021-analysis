# Decision notes

## 2026-10-30: phenotype-to-KASP sample correspondence

The user explicitly confirms E-n corresponds to NYUT_21_n for n = 1 through 53.
Use this user-confirmed numeric-suffix correspondence for the sample mapping.
No independent verification against field records is claimed.

`R/03_sample_mapping.R` reads distinct phenotype entries from the imported CSV
and actual KASP sample metadata from Sheet1!A3:B55. This range excludes headers,
the Indica position row, and later summary rows. The script preserves the exact
KASP ID spelling and Sample Name text and matches numeric suffixes, never row
positions. Repeated phenotype IDs across replications are expected; the mapping
must have unique phenotype IDs, unique KASP IDs, and unique sample names.

Validation on 2026-10-30: 159 phenotype records represent 53 distinct entries;
53 actual KASP sample rows were read. No missing IDs or sample names, duplicate
KASP IDs or sample names, missing expected IDs, or unmatched entries were found.
All 53 entries have a one-to-one mapping. Output: `outputs/sample_id_mapping.csv`.

At the time of this mapping step, the work covered an identifier crosswalk only;
no phenotype models had yet been fitted. Subsequent steps 07–12 fitted provisional
phenotype models and assembled descriptive candidate evidence and QC joins.
No marker–trait association tests have been performed. The correspondence above
remains user-confirmed rather than independently verified from field records.
