# Recalculated scores, full simulations, and manuscript comparison

Date: 21 October 2026. Scope: short- and mid-term pathway analysis in Python and R.

**Completed and independently verified:** both main runs and both historical diagnostic runs exited successfully. There are **16 simulation tables, 33,152 simulated pathway-result rows, and 51,760,000 random-set evaluations**, plus directly calculated hypergeometric results. The recalculated gene scores match the archived values in both languages. **The manuscript's pathway results are not fully reproduced.**

## What was recalculated

Both languages rebuilt 42,575 records from the expression table, gene-ID correspondence, gene lengths, and raw SNP marker blocks. Archived corrected scores were not calculation inputs. The fitted parameters are 0.01 for short-term and 0.47 for mid-term. The exported scores in both languages match **all 42,575 original short values and all 42,575 original middle values exactly as decimal numbers**.

This equality requires explicitly enabled reconstruction hypotheses. Without them, both commands return exit status 2: short has one unresolved infinite score, and middle has two. The other 42,574 short and 42,573 middle scores match the saved values exactly.

| Exceptional calculation | Reconstructed operation | Status |
|---|---|---|
| `Glyma.14G129500`, both phases | Same whole-genome length bin has positive SNP counts 16 and 7; their mean 11.5 rounds to 12. Use the fitted middle parameter and round the corrected probability to seven significant digits. | Numerically consistent; historical operation unconfirmed |
| `Glyma.04G222100`, middle only | Ordinary subtraction produces corrected probability zero. Replace its infinite logarithm with the largest finite middle score, 15.47746851547134, rounded to seven significant digits: 15.47747. | Numerically consistent; historical operation unconfirmed |

The original bin means included zeros. This historical reconstruction is distinct from the general-purpose correction command's newer nonzero-mean/global-fallback rule. The short exception's use of the middle parameter is not established by the author's recollection of averaging. Exact agreement must not be presented as proof of that historical operation. See the [short-score investigation](SHORT_SCORE_RECALCULATION_2026-10-21.md).

## Analysis settings and completed work

The new configurations consume the freshly calculated `score_short` and `score_middle` columns directly, with `score_kind="score"`. They do not apply another logarithm or use the previous WITHOUTINF column-5/6 result files.

| Quantity | Short | Middle |
|---|---:|---:|
| Recalculated records | 42,575 | 42,575 |
| Unique measured genes | 42,567 | 42,567 |
| Primary CTgenes | 44 | 143 |
| CTgenes in the measured universe | 41 | 143 |
| Tested GO terms, measured size 6–1,500 | 2,689 | 2,689 |
| Tested terms containing CTgenes | 98 | 337 |
| Distinct measured pathway sizes | 377 | 377 |

The main workflow averages duplicate scores for the eight extra records. Hypergeometric counts use one common measured universe. Each method applies BH across all 2,689 eligible terms, including terms without CTgenes. SUMSTAT, positive-tail weighted GSEA, and experimental MaxMean each use 10,000 random gene sets per distinct size. Equal-size pathways share a null distribution. Every simulated row records its exceedance count, number of draws, and seed. Seeds are 20261021 and 20261022 for short and middle.

The GSEA implementation ranks genes by decreasing score and reports a positive running-sum enrichment statistic. It is not Broad GSEA's NES/FDR. MaxMean uses the documented global standardization; it does not claim to reproduce every additional transformation in the January 27 exploratory script.

A separate historical SUMSTAT diagnostic retains all 42,575 records, selects full-annotation sizes 5–1,499 with at least one primary CTgene, and uses the old strict-exceedance `b/10000` formula. It tests 100 short and 342 middle pathways, with 95 and 231 distinct record counts. The implementation shares equal-size null distributions; the archived loop sampled separately per pathway, so the marginal null is retained but the joint Monte Carlo dependence and random sequence differ. No original random seed was found.

The main analysis uses `(b+1)/10001`, includes ties, and never reports a zero empirical p-value. Its smallest nominal value is approximately 0.0001. A displayed `1e-16` in an old 10,000-draw table is not evidence of that Monte Carlo resolution.

## Verification and result summary

All 17 staged input hashes match the source manifest. The verifier independently joins GO/CT membership, recalculates duplicate means, checks all saved-score decimals, validates every 10,000-draw count and empirical p-value, and recalculates BH with SciPy independently of the analysis implementation. All result probabilities and statistics are finite.

Python and R agree on memberships, observed statistics, and exact hypergeometric probabilities. The largest observed-statistic difference is 5.01e-12 or less. The largest native Monte Carlo p-value difference is approximately 0.02080; no comparison falls outside the six-standard-error diagnostic screen. The two RNGs need not give identical threshold membership.

All 23 Python tests passed. R numerical, historical-calculation, reference-identity, and end-to-end integration checks passed. The shared-null fixture checks all 120 three-gene subsets of ten genes and confirms identical SUMSTAT, GSEA, MaxMean, finite-null p-values, and FOSCO calculations. Both plotting commands completed, producing 16 figures; representative Python SUMSTAT and R GSEA figures were visually inspected.

Counts below use all 2,689 eligible pathways, before any CT-overlap display filter or secondary biological selection. They are reanalysis outcomes, not the paper's final 13/11 pathway counts.

| Phase | Method | Python BH < 0.05 | R BH < 0.05 | Python BH < 0.005 | R BH < 0.005 |
|---|---|---:|---:|---:|---:|
| Short | Hypergeometric | 9 | 9 | 7 | 7 |
| Short | SUMSTAT | 126 | 126 | 81 | 81 |
| Short | GSEA | 113 | 119 | 74 | 69 |
| Short | MaxMean (experimental) | 126 | 130 | 79 | 82 |
| Middle | Hypergeometric | 139 | 139 | 113 | 113 |
| Middle | SUMSTAT | 247 | 237 | 119 | 113 |
| Middle | GSEA | 124 | 122 | 0 | 0 |
| Middle | MaxMean (experimental) | 228 | 229 | 118 | 113 |

The main runs account for 45,240,000 null-set evaluations; the additional historical SUMSTAT diagnostics account for 6,520,000. Exact hypergeometric tests need no simulation. Runtime records are Python 3.11 / NumPy 2.3.5 / SciPy 1.17.1 and R 4.5.1. R's locale startup notices and experimental-MaxMean notices were nonfatal.

## Manuscript sources and reproducibility limits

The [published article](https://link.springer.com/article/10.1007/s00299-025-03643-2) was checked for the reported cohort, filtering, score preparation, and final pathway selection. Its methods describe fold-change transformation, while the historical scores reconstructed here come from expression p-value columns. The article also describes p-value-based gene statistics in Figure 1. This score-branch discrepancy remains unresolved; the earlier [data-flow audit](HISTORICAL_DATA_FLOW_2026-10-18.md) records the source columns and transformations.

The publisher's two supplementary download URLs returned HTML rather than Office files. No HTML response was accepted as a document. Numerical comparisons therefore use the author's local submission copies in:

`研究法/大豆基因/期刊撰寫/期刊要求的檔案/20250623 plant cell report/Supplementary data/`

- `TableS2_Pathway enrichment analysis.xlsx`: 100 short and 342 middle rows for hypergeometric, SUMSTAT, and GSEA; 100 short and 341 middle MaxMean rows.
- `TableS1.docx`: 13 selected short and 11 selected middle pathways, with three methods each.

Source hashes and individual source ranges are recorded in the verification JSON and comparison TSV. Their identity with the final publisher attachments has **not** been verified. Original files were not modified.

### Hypergeometric discrepancies

The Table S2 values have a specific numerical signature, independently calculated in Python and R:

| Phase | Genuine genes in candidate GO union | Extra distinct pathway-size values | Diagnostic background N | Diagnostic CT total M |
|---|---:|---:|---:|---:|
| Short | 16,860 | 93 | 16,953 | 88 = 2 × 44 |
| Middle | 28,703 | 259 | 28,962 | 286 = 2 × 143 |

These counts follow two operations present in the archived code: counting unique entries from `answer[,3:ncol(answer)]`, which includes the pathway-size column, and taking `length(as.matrix(CTgenes))` on a two-column table. This is a numerical reconstruction of a possible old session state; the saved script also contains column-vector conversions and a different active background expression.

The resulting probabilities match all 100 short Table S2 values within absolute error 5e-10 and all 342 middle values within 5e-8. These intentionally invalid count profiles are saved only as forensic diagnostics. They must not replace the main test's consistent gene universe. The other diagnostic variants change population totals while retaining full-annotation sizes and primary-CT overlaps; they also are not the main coherent-universe test.

Table S1 does not consistently use the same short-term hypergeometric values as Table S2. For example, response to cold is 3.78e-10 in S1 and 0.01219794 in S2. Using S1's printed pathway/CT counts, most selected short values match a background of 42,567 with 41 CTgenes; the remaining two match the older container-count profile. The machine-readable report lists the profile matches individually.

### SUMSTAT and GSEA discrepancies

Even the historical SUMSTAT replay differs materially from Table S2: 36 of 100 short comparisons fall outside a broad six-standard-error Monte Carlo screen in each language; middle has 112 of 342 in Python and 114 in R. This is a diagnostic screen, not an equivalence test. The differences cannot all be explained by the native random generators or by averaging duplicates in the main workflow. A consistent archived execution state tying the exact input branch to the manuscript table has not been established.

There are also discrepancies between the submission tables themselves. Short-term response to water has SUMSTAT p = 0.0171 in S1 and 0.7543 in S2. Regulation of gibberellin biosynthetic process has GSEA p < 1e-16 in S1 and p = 0.0324 in S2. These are not ordinary printed rounding differences.

All 100 short and 342 middle S2 SUMSTAT p-values match the saved `篩選後/SUMSTAT/SUMSTAT_ct_*_pavlue.txt` files within 1e-12. The discrepancy therefore precedes this table extraction. The adjacent `篩選後/SUMSTAT/程式.txt` uses another branch: phase-specific subsets of `outside_data`, a different infinite-score replacement, and, for middle, different objects for the observed and null scores. These saved scripts depend on session variables, so this evidence does not establish a complete executable provenance chain for S2.

Archived GSEA versions sort by gene identifier, use differing observed/null inputs, or depend on previously saved ES values. The cleaned score-ranked GSEA is therefore an explicitly different analysis. Its successful execution does not validate those older GSEA p-values. The full TSV compares every available pathway rather than selecting only agreements.

### Pathway eligibility and final selection

The main workflow can compare 97 of the 100 short and 335 of the 342 middle S2 pathways. The remaining terms fail its measured-size bounds. Of the final 13 short pathways, three are excluded: response to water (`GO:0009415`), regulation of gibberellin biosynthetic process (`GO:0010371`), and multicellular organism growth (`GO:0035264`); their full annotations contain only five genes. All 11 selected middle pathways are eligible.

The main BH family contains 2,689 terms, whereas the historical candidate family contains 100/342. Changing this family changes adjusted probabilities. Final published pathway selection is not forced to 13/11, and the reported count/ratio/top-fraction selection rules are not silently changed to obtain those totals.

## Commands and output locations

The [README](../README.md#recalculate-both-phases-and-run-full-simulations) gives both recalculation, analysis, plotting, and historical diagnostic commands. All analysis outputs are under `results/recalculated_2026-10-21/`:

| Directory | Contents |
|---|---|
| `python_scores`, `r_scores` | Fresh two-phase scores, fitted grids, donor/exception audits, and provenance |
| `python_strict`, `r_strict` | Baseline calculations without exceptional repair hypotheses; expected exit status 2 |
| `python_full`, `r_full` | Main method tables, unique-gene scores, membership audits, settings, and plots |
| `python_historical`, `r_historical` | Historical SUMSTAT simulations and explicitly labeled hypergeometric variants |
| `publisher_sources` | Failed Office-download status, including the actual HTML content type |

For the independent baseline comparison, also run both recalculation commands without the hypothesis flag. Each writes its audit outputs and deliberately exits with status 2 because unresolved infinite scores remain:

```sh
python scripts/recalculate_phase_scores.py "D:/Bioinformatics project" results/recalculated_2026-10-21/python_strict
Rscript --vanilla scripts/recalculate_phase_scores.R "D:/Bioinformatics project" results/recalculated_2026-10-21/r_strict
```

The comparison extraction requires `openpyxl`; runtime analysis itself does not. Replace `SUPPLEMENT_DIRECTORY` below with the local manuscript directory listed above:

```sh
python -m pip install openpyxl
python scripts/extract_manuscript_enrichment.py "SUPPLEMENT_DIRECTORY/TableS2_Pathway enrichment analysis.xlsx" results/recalculated_2026-10-21/manuscript_reference.json --selected-docx "SUPPLEMENT_DIRECTORY/TableS1.docx"
python scripts/verify_recalculated_runs.py --source-root "D:/Bioinformatics project"
```

The verifier independently reads archived reference scores from `--source-root`. This run also records the observed zero process-exit codes and file timestamps in `execution_record.json`; a future CLI run can validate saved outputs without that optional execution record.

The full numeric comparison is [recalculated_pathway_comparison_2026-10-21.tsv](recalculated_pathway_comparison_2026-10-21.tsv). Verification, hashes, method summaries, and cross-language differences are in [recalculated_full_verification_2026-10-21.json](recalculated_full_verification_2026-10-21.json). Generated result folders remain local and are ignored by Git.
