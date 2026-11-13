# Actual input audit — 2026-11-09

These checks were executed with Python's gzip/csv standard-library readers, independently of the unexecuted R script. All supplied files decompressed and parsed successfully.

- Counts: 39,376 rows; six numeric sample columns.
- Missing, non-finite, negative or noninteger counts: none found.
- Duplicate count GeneIDs: none found.
- Annotation: 39,376 unique GeneIDs in exactly the same order as counts; all count rows matched.
- FPKM and TPM: 39,376 rows each; same six sample column names. Not used as input for statistical testing.
- Series matrix: 12 sample records; the six RNA-seq libraries match the count columns. The other six are labelled miRNA-seq and are excluded.
- The series matrix has an empty expression table: it is used for metadata only.
- Donor IDs, pairing and replicate independence are not established by this metadata file.

| Sample | Day | Assigned counts | Genes with count > 0 |
|---|---|---:|---:|
| GSM5395968 | 0 | 17,513,897 | 25,357 |
| GSM5395969 | 0 | 27,407,175 | 23,026 |
| GSM5395970 | 0 | 23,791,049 | 23,566 |
| GSM5395971 | 21 | 26,283,576 | 24,441 |
| GSM5395972 | 21 | 25,443,410 | 25,660 |
| GSM5395973 | 21 | 18,188,936 | 24,356 |

The prespecified descriptive filter (CPM >= 1 in at least three libraries) retains **15,562 genes**. This was computed from actual counts, not assumed. It is not the number of differentially expressed genes.

Unequal library sizes are observed. This alone is not evidence of failed samples. No sample was removed. Read quality, mapping rates and contamination cannot be inferred from this count-table audit.

No PCA coordinates, plots, p-values, differential-expression results or enrichment statistics have yet been generated with R. Run `run_qc.R` in RStudio for the next outputs.
