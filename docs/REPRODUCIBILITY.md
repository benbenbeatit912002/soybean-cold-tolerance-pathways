# Reproducibility status

This repository is a research companion and documented reanalysis. Running without errors is different from reproducing the published results.

## Score reconstruction

The portable workflow rebuilds 42,575 records from four packaged upstream sources. It does not read expected scores during calculation. Duplicate IDs remain distinct during reconstruction; pathway analysis subsequently averages scores for the 42,567 unique genes.

The strict archived calculation and inferred export precision reproduce 42,574 short-term and 42,573 mid-term saved scores. Remaining nonfinite values stop the pipeline with exit status 2. The explicit `--use-reconstruction-hypotheses` option enables:

1. **The remaining zero-SNP record, Glyma.14G129500.** Positive whole-genome donors in its length bin have counts 16 and 7; their mean 11.5 rounds to 12. The short repair uses the fitted middle-phase lambda 0.47 and seven-significant-digit probability rounding. Averaging was reported as the historical repair method; lambda reuse and the exact operation remain unconfirmed.
2. **The middle-phase cancellation case, Glyma.04G222100.** A probability rounded to zero under historical arithmetic is represented by a score calculated from the largest finite score and rounded to seven significant digits. Its historical use is unconfirmed.

With these hypotheses, both languages calculate all 42,575 archived scores in each phase exactly at the saved decimal precision. This numerical agreement does not establish that the inferred repairs were the actual historical procedure. Read the [full score and pathway rerun report](../reports/RECALCULATED_ANALYSIS_2026-10-21.md).

The historical imputation rule includes zeros in the group mean and has no global fallback. The separate `scripts/correct_scores.*` utilities implement the alternative nonzero group mean/global nonzero mean rule with stable arithmetic. These are distinct procedures. The portable main workflow follows the reconstructed historical branch to make that distinction inspectable.

## Input lineage and Yamasaki workbooks

The included expression source is a byte-preserved copy of the converted Peixiu TXT table. Its columns 5 and 8 supply short- and mid-term expression p-values. The original external-import script reads that TXT, not an Excel workbook.

Yamasaki workbooks were inspected during provenance work, but the complete workbook-to-TXT conversion chain is not established. The portable package starts from the supplied TXT and does not claim to regenerate it from Excel. Fold-change columns exist in other branches; they are not substituted into this p-value branch.

## Pathway analysis

The main workflow fixes these choices explicitly:

- Short and middle are separate; CT prior scores define stored annotations and are not added to expression scores.
- The background contains 42,567 measured unique genes; eight extra duplicate records are averaged.
- GO eligibility uses measured sizes 6 through 1,500, producing 2,689 tested terms per phase.
- SUMSTAT and GSEA use 10,000 random sets per distinct measured size in full mode. MaxMean is optional.
- A simulation p-value is `(exceedances + 1)/(draws + 1)`, counting null statistics at least as large as observed.
- BH adjustment includes the complete eligible family, separately by phase and method.

The paper reports 42,576 in its workflow, uses different candidate term counts, and describes a fold-change transformation. Archived scripts and result tables also contain inconsistent background/counting and GSEA variants. Selected local manuscript tables are not entirely mutually consistent. Consequently, matching all published pathway p-values is **not verified**.

The [October 30 membership audit](../reports/CT_MEMBERSHIP_LINEAGE_2026-10-30.md) resolves the two selected short CT-count differences through the 44-to-41 measured foreground. Python and R reproduce the 100 archived January short hypergeometric rows at saved precision; R also reproduces the 342 middle rows. Stable direct-tail values agree across languages, while a few middle legacy subtraction values differ at machine precision in Python. Two final short Excel rows still use an older calculation version. These findings refine the historical audit without changing the main statistical workflow or confirming the candidate selection rules.

The [paper concordance audit](../reports/PAPER_CONCORDANCE_2026-10-18.md) and [full comparison](../reports/RECALCULATED_ANALYSIS_2026-10-21.md) record the evidence. Local manuscript files were compared; final publisher supplement identity was not verified. Historical mistake-replay configurations are diagnostic tools, not the recommended main analysis.

## What tests establish

The [crosstalk extension](CROSSTALK.md) uses the author's historical Excel selection, verified against the saved 13/11 lists and local Table S1. The paper documents the selection considerations. A configurable selector reproduces the IDs from full historical candidate tables, but its exact rule remains unconfirmed. Recomputed Jaccard coefficients match the archived matrices at printed precision. The recovered dated script's cutoffs and nine-decimal comparison reproduce 61 edges; the earlier fraction interpretation gives 67. The middle 121/122 difference is traced to Glyma.09G149200 and a saved 14-pathway gene export versus the final 11-pathway list. See the [source trace](../reports/CROSSTALK_TRACE_2026-10-24.md).

Tests check formulas, input identity, duplicate handling, imputation, shared-null cross-language agreement, executable entry points, and independent verification of outputs. Portable package tests start from a source ZIP with no saved results, in a new directory containing spaces.

Python and R use different random generators. Native simulation p-values can differ; they are not required to be identical. For every run, retain the saved settings and language environment. Local Windows validation does not establish that the configured GitHub Actions jobs have already passed on remote runners.
