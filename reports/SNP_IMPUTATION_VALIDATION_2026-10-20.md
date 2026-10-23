# SNP imputation and correction validation — 20 October 2026

**Follow-up, 21 October 2026:** the author requires the original short scores as the paper reference. The [reference-preserving rerun](PAPER_SHORT_REFERENCE_VALIDATION_2026-10-21.md) now retains all saved short scores exactly and keeps the recalculated values in separate columns. The results below remain a record of the earlier reanalysis; its changed short scores are not the authoritative paper scores.

**Result: PASS.** Python and R implement the nonzero-group-mean/global-fallback replacement rule and completed a real-data FOSCO rerun for short and middle phases. All 42,575 records were compared between languages, against the saved input, and against an independent calculation. The single remaining zero SNP count became **11**. The largest absolute difference between Python and R's final `-log10` scores was **4.98e-14**.

This validates imputation and correction from the saved historical input. It does not establish exact reproduction of the paper or rerun the full pathway simulations. Machine-readable evidence is in [snp_imputation_verification_2026-10-20.json](snp_imputation_verification_2026-10-20.json).

## Was this the original method?

Not exactly. The original import script used a rounded group mean, but **included zeros in that mean** and **had no global fallback**. Its interval loop also ran beyond valid breakpoints. See [`legacy/#外部資料匯入與校正繪圖.txt`](../legacy/%23外部資料匯入與校正繪圖.txt), lines 61–68.

| Decision | Original script | Updated Python and R |
|---|---|---|
| Grouping variable | Gene length, 2,000-bp intervals | Same grouping, with complete bounds |
| Donors in a group | All stored counts, including zero | Original positive input counts only |
| No positive donors in the group | Group mean can remain zero | Overall mean of original positive counts |
| Timing of means | Mean calculated once for each group | Group and global means use original positive counts |
| Rounding | R `round()` | Nearest integer, ties to even, in both languages |
| Record of each replacement | No standalone audit | Original count, new count, bin, donors, mean, and source |

Groups are `[1,2001)`, `[2001,4001)`, and so on, without an upper length limit. Missing, negative, or fractional SNP counts are rejected. An all-zero input fails because there is no positive donor population.

The new rule is implemented in [`python/ctpathway.py`](../python/ctpathway.py) and [`R/core.R`](../R/core.R), and called by both `scripts/correct_scores` entry points before fitting lambda.

## Input and replacement

The source is `研究法/大豆基因/數據/外部基因校正/gene_p_len.txt` under the original project root. The staging script preserves its raw short/middle values, lengths, counts, order, and all records, while naming columns explicitly.

| Check | Result |
|---|---:|
| Input records | 42,575 |
| Unique canonical gene IDs | 42,567 |
| Extra duplicate records preserved | 8 |
| Zero SNP records | 1 |
| Original positive records used as global donors | 42,574 |
| Global nonzero mean | 10.731667214731996 |
| Rounded replacement | 11 |

The affected gene is **`Glyma.14G129500`**, at physical source line **29,915** (data row 29,914). Its length is **52,247 bp**, placing it in `[52001,54001)`. That group contains no positive donor count, so the global fallback applies. The raw SNP-set block contains this gene followed immediately by `END`, supporting a count of zero assigned markers in that file; it does not establish that the gene has no biological variation.

**Scope of the donor population:** `gene_p_len.txt` already contains earlier historical imputations. This run changes its remaining zero and preserves its stored positive counts. It does not recompute all original zero-marker cases from raw SNP records.

Both commands explicitly use `--duplicate-policy preserve` for this historical comparison. Thus donor means and lambda fits count input records, including the eight extra rows. The default rejects duplicate genes. Preserving these records for correction does not settle their aggregation policy for future gene-level pathway analysis.

## Numerical verification

| Check | Short | Middle |
|---|---:|---:|
| Selected lambda | 0.01 | 0.47 |
| Maximum Python/R difference in adjusted values | 6.11e-16 | 6.67e-16 |
| Maximum Python/R difference in `-log10` scores | 4.98e-14 | 4.98e-14 |
| Maximum difference from independent stable formula, Python | 0 | 0 |
| Maximum old/new adjusted-value difference among the 42,574 originally positive-SNP records | 7.22e-16 | 7.22e-16 |
| Old/new adjusted-value differences greater than 1e-12 | 1 | 1 |

Both languages chose the same lambdas as the recovered historical correction. All 100 grid points per phase were compared. All **85,150 phase-specific corrected values are positive** and all final scores are finite.

The independent reference uses Python's standard-library `math.expm1` and `math.log1p`, separately from the NumPy implementation. The verifier also checks source-row identity, gene order, lengths, raw scores, unchanged positive SNP counts, independent donor means, and output hashes.

Completed checks:

- All **14 Python tests** passed, including the synthetic command-line analysis, correction, and plotting workflow.
- R's numerical and synthetic integration checks passed.
- Imputation tests cover group/nonzero means, empty-group fallback, frozen donors, row-order independence, bin boundaries, lengths above 100,000 bp, ties-to-even rounding, and invalid inputs.
- The cross-language reference passed for SUMSTAT, GSEA, experimental MaxMean, identical-null p-values over all **120** three-gene subsets of a ten-gene universe, and fixed/fitted FOSCO.
- Both real-data correction commands exited successfully and generated short/middle diagnostics.
- The new Python short-phase and R middle-phase diagnostic plots were visually inspected; labels and panels render without clipping.

## Comparison with the historical values

For `Glyma.14G129500`, the old full-precision correction file stores SNP **12**; the new rule produces **11**:

| Phase | Old adjusted value | New adjusted value | `-log10` of old stored value | New score |
|---|---:|---:|---:|---:|
| Short | 0.9947477 | 0.812162816978 | 0.002287056332 | 0.090356897573 |
| Middle | 0.9008147 | 0.891192146373 | 0.045364535374 | 0.050028649372 |

The large short-score difference cannot be explained solely by replacing 12 with 11. Applying the same stable formula to the verified raw short value, lambda 0.01, and SNP 12 produces score **0.090210820639**, not the old stored score **0.002287056332**. The isolated 12-to-11 effect is only **+0.000146076934** in the short score. The historical exceptional value remains unexplained.

At lambda 0.47, the isolated 12-to-11 effect is **+0.004664095825** in the middle score. This is consistent with the scale of the observed middle change; the old stored value has limited decimal precision.

One separate numerical issue is visible for **`Glyma.04G222100`**, whose SNP count remains 9. Its raw middle value is `5.21e-17`. Stable arithmetic returns adjusted value **1.4632943736267197e-16** and score **15.834668297264**. The historical correction stores zero, while the historical logged file separately stores **15.47747**. The improvement from zero to a positive corrected value comes from the existing stable formula, not from the new SNP replacement rule. The old finite logged-score substitution remains unexplained. A tiny absolute difference in a probability can therefore matter greatly after `-log10`.

Full exceptional-row details are saved locally in `results/snp_imputation_2026-10-20/historical_gene_comparison.tsv`, and included in the JSON evidence.

## Relationship to pathway results

Imputation preserves all gene IDs and memberships. Hypergeometric results are therefore unchanged **when the same universe, CT membership, pathway rules, and testing family are used**. The imputed gene belongs to `GO:0009507`, `GO:0009630`, `GO:0009793`, and `GO:0048046`. This does not imply that score-based changes are confined to these four terms: GSEA ranking, null distributions, and MaxMean standardization can depend on the whole score universe.

The prior 10,000-draw configurations analyze WITHOUTINF columns 5 and 6 directly. They do not read `gene_p_len.txt` or these correction outputs. Their existing results are separate execution records and must not be labelled a rerun with the new SNP policy. This validation does **not** claim new SUMSTAT/GSEA/MaxMean pathway p-values from the corrected historical branch.

The distinction between the historical p-value input and the paper's described fold-change input, duplicate handling for pathway sampling, and final paper selection settings still requires reconciliation; see the [paper concordance audit](PAPER_CONCORDANCE_2026-10-18.md).

## Reproduce this correction comparison

Run from the repository root. The historical source files must exist under the supplied source root; they are not staged by the ordinary WITHOUTINF preparation command.

```sh
python scripts/prepare_historical_correction.py "D:/Bioinformatics project" results/snp_imputation_2026-10-20/input.tsv
python scripts/correct_scores.py results/snp_imputation_2026-10-20/input.tsv results/snp_imputation_2026-10-20/python --duplicate-policy preserve
Rscript --vanilla scripts/correct_scores.R results/snp_imputation_2026-10-20/input.tsv results/snp_imputation_2026-10-20/r --duplicate-policy preserve
python scripts/verify_snp_correction.py "D:/Bioinformatics project" results/snp_imputation_2026-10-20
```

```sh
python -m unittest discover -s tests -p "test_*.py"
Rscript --vanilla tests/test_core.R
Rscript --vanilla tests/test_integration.R
python tests/check_cross_language.py
```

On the author's Windows installation, invoke R with `& 'C:/Program Files/R/R-4.5.1/bin/Rscript.exe'` in PowerShell if it is not on PATH. The cross-language checker accepts `--rscript "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"`. Validation used Python 3.11 and R 4.5.1; output directories contain the recorded package/runtime versions.

The source files and historical outputs were preserved. Current changes and generated results are local; this report does not assert that they have been pushed to GitHub or added to the older review ZIP.
