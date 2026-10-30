# CT membership lineage and historical hypergeometric replay

The two short-term CTgene count differences are explained by a change from the complete 44-gene primary list to the 41 genes retained in the historical expression universe. Both implementations independently recover the archived counts from gene IDs. No gene or pathway probability was inserted to force agreement.

## The two selected pathways

| GO ID | Pathway | Primary CT count | Measured CT count | Final Excel count | Excluded gene |
|---|---|---:|---:|---:|---|
| GO:0045893 | Positive regulation of transcription, DNA-dependent | 15 | 14 | 14 | Glyma.01G216200 |
| GO:0050832 | Defense response to fungus | 12 | 11 | 11 | Glyma.01G216200 |

The complete historical candidate workbook uses the primary counts 15 and 12 in all three main methods. The original final Excel uses 14 and 11. The archived `short_CT_G_score.txt` contains exactly the intersection of the primary short CT list with the 42,567 unique mapped expression genes. The original hypergeometric script loads that scored CT file at lines 11-13, computes overlap at lines 39-40 and 72, and exports the overlap at line 83.

The current main Python and R analyses already use measured CT membership and produce 14 and 11. The candidate Excel selector preserves its archived input counts. These are distinct outputs with distinct input definitions; the selector is not a fresh enrichment run.

All 24 final selected rows (13 short and 11 middle) have the same CT counts as a fresh intersection with measured CT membership. The original 4,896 GO membership sets also agree exactly with the packaged GMT; no annotation change is needed to explain these two rows.

## Why three primary short CTgenes are absent

| Gene | Original alias | Evidence in the supplied inputs |
|---|---|---|
| Glyma.01G216200 | Glyma01g42510 | A correspondence exists at row 1977, but this alias has no record in the expression table. |
| Glyma.12G182200 | Glyma12g30710 | A correspondence exists at row 30040, but this alias has no record in the expression table. |
| Glyma.07G159000 | Glyma07g19702 | Expression row 15814 matches correspondence rows 16365 and 16366. Both point to the same modern gene, but the original code requires exactly one matching correspondence row. |

Rows are physical file lines, including the header. The original import script's `length(loc)==1` condition is at line 46. Thus the third exclusion comes from the row-count rule, even though the modern target is unambiguous. This finding does not establish that an absent expression record means zero biological expression. The audit leaves the historical import policy unchanged.

Only Glyma.01G216200 belongs to the two selected short pathways with differing counts. The other two excluded genes do not explain those rows.

## Independent replay and its limits

Both audit commands calculate the old candidate family using full annotation sizes 5 through 1,499 and at least one measured CTgene. They derive all counts from IDs, use N = 42,567 unique measured genes and M = 41 short / 143 middle CTgenes, and apply BH across the 100 / 342 historical candidates. Expected table values are read only by the verification test after calculation.

| Check against original January 09 tables | Python | R |
|---|---|---|
| Short: all 100 CT counts, nominal p and BH values | Match saved precision | Match saved precision |
| Middle: all 342 CT counts | Exact | Exact |
| Middle: nominal p and BH using the old subtraction expression | Three nominal cells and four BH cells differ near machine epsilon | All 342 rows match saved precision |
| Stable direct-tail probabilities and BH across both languages | Agree within relative tolerance 2e-11 | Agree within relative tolerance 2e-11 |

The small Python middle differences arise from subtracting a CDF near one. They affect GO:0019288, GO:0019216 and GO:0015994 in nominal probabilities, plus GO:0031408 through BH adjustment. Both commands retain stable direct-tail results and the legacy subtraction results separately. No arbitrary zero replacement is applied.

This advances reproduction of the archived January tables, not of every final manuscript value. In the final Excel, 11 of 13 short hypergeometric rows follow the January calculation. Photoprotection (GO:0010117) and carbon utilization (GO:0015976) retain values consistent with the older population-count profile described in [the recalculated analysis report](RECALCULATED_ANALYSIS_2026-10-21.md). Their final Excel p-values are 0.065440935 and 0.060561021, versus 0.012451074 and 0.011498693 under the January convention. The final table therefore mixes calculation versions.

The main reanalysis continues to restrict pathway membership to the measured universe and uses its own eligible testing family. Historical annotation-size p-values must not silently replace the main results. The candidate secondary-selection rule remains unconfirmed.

## Remaining mid-term candidate count inconsistency

The complete candidate workbook has 17 mid-term pathways whose SUMSTAT and GSEA CT counts are each one lower than the primary/measured intersection. The hypergeometric counts agree. All 143 primary middle CTgenes are present in the mapped expression universe, so the short-term exclusion mechanism cannot explain these 34 method rows.

Deleting any one of six CTgenes produces the same 17-pathway count pattern: Glyma.02G005600, Glyma.03G250600, Glyma.10G006600, Glyma.10G165800, Glyma.10G180800 or Glyma.20G209700. This pattern does not identify which gene, if any, was omitted in the old session. The saved scored middle CT file contains all 143 genes. No specific omitted gene or operation is asserted. This is separate from the previously resolved 14-versus-11 pathway version difference behind the 122-versus-121 network export.

## Run the audit

After a normal quick or full pipeline run, supply its freshly recalculated gene table. Only `gene_id` is used; all score columns are ignored.

```sh
python scripts/audit_historical_membership.py results/python_quick/scores/recalculated_scores.tsv results/membership_python
Rscript --vanilla scripts/audit_historical_membership.R results/r_quick/scores/recalculated_scores.tsv results/membership_r
```

Use a new output directory each time. Each command exports per-pathway primary and measured CT counts, excluded IDs, the historical hypergeometric replay, input checksums and a summary. They work from other directories when given absolute input and output paths. No external network service is used.

```sh
python -m unittest discover -s tests -p test_historical_membership.py
```

The automated test uses the gene-ID column of the packaged historical input for a fast deterministic check. It compares expected IDs and saved probabilities only after calculation; the archived gene scores are not used.

## Evidence and scope

The source audit read the original research files without changing them. It checked all 4,896 annotation sets, all 24 final selected CT counts, and all 43,024 CT-count entries across the existing Python/R full-run method tables. SUMSTAT, GSEA and MaxMean simulations were not rerun for this membership-only change.

Original paths below are relative to the research archive:

- `CTgenes 含分數/SHORT.txt` and `數據/CT_genes/short_CT_G_score.txt`.
- `數據/GO_terms_4899.txt`.
- `程式/超幾何.txt` and `程式/#外部資料匯入與校正繪圖.txt`.
- `篩選後/篩選後大豆ct short0109.txt` and `篩選後/篩選後大豆ct middle0109.txt`.
- `期刊撰寫/table/Table 1. The results of the pathway enrichment analysis with secondary selection..xlsx`, Short-term!A8:O20 and Mid-term!B10:J20.

Validation snapshots and source hashes are in [historical_membership](../data/reference/historical_membership/README.md). The machine-readable findings are in [ct_membership_verification_2026-10-30.json](ct_membership_verification_2026-10-30.json). Local June submission tables remain distinct from publisher attachments whose identity has not been verified.
