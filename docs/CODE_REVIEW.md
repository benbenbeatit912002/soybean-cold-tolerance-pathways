# Code and data review — 15 October 2026

**Recalculation update, 21 October 2026:** both historical calculators rebuild inputs from four upstream files. Saved scores serve only as independent verification references; the earlier copying option was removed. The archived calculation and inferred export rules match 42,574 saved short scores exactly. A separately enabled repair hypothesis computes the remaining value from the whole-genome bin's positive mean, the fitted middle lambda, and seven-significant-digit rounding. It matches all 42,575 scores in both languages, but the original repair operation is still unconfirmed. See the [recalculation report](../reports/SHORT_SCORE_RECALCULATION_2026-10-21.md).

**Evidence update, 18 October 2026:** see the [historical data-flow audit](../reports/HISTORICAL_DATA_FLOW_2026-10-18.md) and [script dependency inventory](../reports/SCRIPT_DATA_DEPENDENCIES_2026-10-18.md). The original mapping and both saved FOSCO transformations have now been numerically reconstructed, apart from a zero-SNP exception. Most previously noted logged-score discrepancies are explained by comparing against a rounded intermediate; the high-precision corrected source reproduces the finite logged scores within `5e-8`, with two documented exceptional genes. Earlier unresolved findings below should be read with this newer evidence.

**Implementation update, 20 October 2026:** the alternative nonzero-group-mean/global-fallback SNP policy is implemented in Python and R. A [complete correction rerun](../reports/SNP_IMPUTATION_VALIDATION_2026-10-20.md) from the saved historical input passed comparison across languages and against independent calculations. The new rule resolves the remaining zero for this reanalysis; it does not explain how the original saved exception was produced.

## Assessment

The supplied material can be organized as a research repository, but the original scripts cannot yet support an exact end-to-end reproduction claim. The audit found input-version mismatches, unresolved gene mappings, parse failures, and statistical implementation differences. Original files were preserved. A checked reanalysis workflow and explicit data lineage are included alongside them.

All 11 files in `研究法/大豆基因/程式` were read and parsed without execution. Additional checks covered the primary CTgene lists, the two specified WITHOUTINF inputs, GO annotation, correspondence data, stored correction intermediates, and representative upstream/server scripts. Existing published figures and supplementary result tables were not exhaustively reproduced or independently validated.

Priority meanings: **P1** can change results or prevent execution; **P2** affects reproducibility, interpretation, or figure accuracy. Line references refer to the unmodified files in `legacy/` unless another source is named.

## P1 findings

### 1. The current statistics scripts use a different external-data branch

`超幾何.txt:6–15`, `sumstat.txt:2–4`, and `maxmean0127版本.txt:5–7` load `Yamasaki_Go校正後1213_G_log10.txt`, not the two specified WITHOUTINF files. The corrected file has 42,575 rows and 42,567 unique IDs; the WITHOUTINF files have 40,061 and 40,332 raw rows. `short_CT_G_score.txt` has 41 genes rather than the primary list's 44.

**Consequence:** changing the README alone would misstate the analyzed inputs. **Action:** stage the author's named primary and external inputs, retain hashes, and label new results as a reanalysis. Historical output attribution remains unresolved.

### 2. The short external file mixes genome-ID versions and both files contain blank IDs

Although named `2.0`, the short file has 1,396 older-format IDs. The supplied correspondence table resolves 1,388 to unique v2 targets; eight have multiple targets. There are five blank-ID rows in short and two in middle. These are distinct from numeric `Inf` checks.

**Consequence:** case normalization alone misses valid genes; blindly counting rows can count missing identifiers or the same gene more than once. **Action implemented:** remap unique targets, exclude and record ambiguous/blank identifiers, retain unplaced v2 genes, and export every row's decision. Resolved unique universes are 39,958 and 40,235 genes. See `reports/import_summary.json` and the two import-decision TSVs.

### 3. Repeated gene IDs have different scores

After mapping and removing invalid identifiers, the external files contain 90 and 95 extra rows beyond their unique gene counts. The stored correction branch also contains eight duplicated genes with differing values. Thus the duplicates are not generally identical rows that can be discarded without a scientific choice.

**Consequence:** row-level sampling overweights repeated genes; `%in%` can bring multiple scores into a supposedly gene-level statistic. **Action implemented:** membership tests use unique IDs. Score analyses reject duplicates by default. Diagnostic configurations explicitly average the supplied values per gene; this is a sensitivity assumption, not a reconstruction of the author's historical policy. Resolve whether these rows represent alternative source mappings, transcripts, or independent measurements before a final rerun.

### 4. Two scripts fail parsing; several others rely on a previous interactive session

R 4.5.1 parsing fails at `超幾何.txt:126` (an uncommented separator) and `比較分數與無母分析程式.txt:132` (a dangling comma and parenthesis). See `reports/legacy_parse.tsv`.

Even after removing the separator, `超幾何.txt:13` converts `CTgenes_short` into a vector, then line 39 uses `CTgenes_short[,1]`, causing an incorrect-dimensions error. `CTgenes_middle` has the same issue at lines 16 and 50. Line 113 refers to undefined `answer4`. The file contains a short-term probability loop but no completed equivalent middle-term loop.

`sumstat.txt` assumes `answer`, `answer2`, `k1`, `k2`, and selection indices already exist. MaxMean assumes matching overlap objects. `gsea全部.txt:91` and 118 contain bare symbols, which parse but will fail unless such objects happen to exist; line 49 reads `datfile` instead of the newly assigned output filename. Parsing success is therefore not execution success.

**Action implemented:** explicit entry points and function arguments; no reliance on an RStudio workspace. Legacy scripts are not executed automatically.

### 5. Hypergeometric counts do not consistently describe one universe

`超幾何.txt:63–73` takes L from the corrected-score IDs but uses each GO term's original annotation size for s and the full CT list for M. Those sets are not explicitly intersected with L. Pathway selection occurs earlier using a CT list from another file.

**Consequence:** the population, successes, draws, and overlap can refer to different sets of genes. **Action implemented:** intersect every set with the same resolved measured universe and recompute all counts. Missing CTgenes are exported. Use the direct upper tail `phyper(..., lower.tail=FALSE)` instead of `1-phyper(...)` to avoid numerical cancellation. The comment at line 150 recommends the less stable direction.

### 6. GSEA ranks gene names, and some permutation code compares different objects

`gsea_菲斯福氣.txt:50`, 99, 160, and 208 and `gsea全部.txt:9` and 64 order column 1, which contains IDs. The representative server `研究法/short_coding.R:24` also orders IDs. GSEA should rank the intended numerical metric.

In `gsea_菲斯福氣.txt:47–48`, the observed statistic uses the CT-only middle list; lines 95–99 switch the null to all genes. Lines 118–121 sample numerical score values and then use them as row indices. Lines 129–132 still traverse the unshuffled object. The later test section has related inconsistencies. `gsea全部.txt:97–98` repeats the score-as-index problem.

**Action implemented:** one score universe for observed and null statistics, numerical sorting, and an independently checked positive-tail running sum. This is a new GSEA-style implementation with documented limits, not a relabelled copy of the historical output.

### 7. Permutation p-values omit ties, can be zero, and are not reproducible

`sumstat.txt:48–56`, 82–90 and `maxmean0127版本.txt:54–55`, 116–117 count strict exceedances and accumulate `1/B` without a seed or finite-sampling correction. GSEA files contain several different test iteration counts, including 2, 10, and 100 in exploratory code. Server files save exceedance counts whose pooling and denominator require tracing.

**Action implemented:** seeded random sampling, include ties, and `(b+1)/(B+1)`. Keep nominal and BH values separately and record B. Previously reported zero values need rerunning, not replacement with an arbitrary tiny number. See [Phipson and Smyth](https://pubmed.ncbi.nlm.nih.gov/21044043/).

### 8. The correction input and corrected phase are not consistently identified

`#外部資料匯入與校正繪圖.txt:50–51` reads columns 5 and 8 of the transformed upstream table, whose headers identify them as p-value fields. The [paper](https://doi.org/10.1007/s00299-025-03643-2) describes a fold-change-derived transformation. This difference cannot be settled from variable names alone.

At lines 150–170, an object named `p_middle` is created, but only column 2 (short term) is adjusted and written. This exact script does not apply a separate fit to column 3. The plotting script `校正前後的繪圖.txt:19–33` also changes short only, while later displaying middle as corrected. Its objective at line 29 is gene length; the import script's objective at line 165 is SNP count. The upstream `研究法/外部資料處理.txt:103` uses Pearson instead of Spearman.

**Action implemented:** separate per-phase correction functions, explicit required inputs, stable arithmetic, and recorded objectives. The 20 October rerun uses the actual saved historical p-value columns and reproduces the ordinary corrected records; the choice between this branch and the paper's fold-change-derived branch still needs reconciliation. The supplied general external scores are not silently corrected.

### 9. Zero-SNP imputation has uncovered boundaries and no explicit failure policy

`#外部資料匯入與校正繪圖.txt:63–68` constructs breakpoints ending at 98,001 and iterates past valid interval indices. The last open-ended range is not represented. Means include zero counts, and rounding can leave zeros. A zero M makes the correction return zero for ordinary inputs, leading to an infinite negative log.

The saved `gene_p_len.txt` contains one zero SNP count; the corrected saved file has none. The step that resolves this difference is not fully recorded in the reviewed script.

**Action implemented, 20 October 2026:** zero counts use the same 2,000-bp group's nonzero mean, falling back to the overall nonzero mean. Bins extend to all positive lengths, means use positive input records before any imputation, and rounding uses nearest integer with ties to even. Missing/negative/fractional counts and all-zero inputs are rejected. Original counts and every replacement are recorded. The saved input's sole zero, `Glyma.14G129500`, has no nonzero donors in `[52001,54001)` and becomes 11 from the global mean `10.731667214731996`. Its old saved SNP value 12 and exceptional short score remain historical differences; the new rule should not be described as the original rule.

### 10. January 27 MaxMean adds an unverified within-set transformation

`maxmean0127版本.txt:15–36` computes global z-scores and positive/negative means. Lines 39 and 53 then center and divide the MaxMean statistic by each pathway's or sample's own mean and SD, with corresponding steps at lines 101 and 115. This changes the statistic and can divide by zero; it is not established as the intended restandardization procedure by the code or comments. `maxmean原版.txt:25` and 33 instead calculate `(sum - mean)/sd`, a different statistic again.

**Action implemented:** keep both originals, provide a clearly experimental conventional MaxMean statistic, and do not run it by default. An unsatisfactory result alone does not establish an implementation error; the extra transformation warrants methodological review independently of its outcome.

## P2 findings and interpretation limits

| Finding | Evidence | Action |
|---|---|---|
| Pathway bounds differ from the stated methods | `超幾何.txt:28`: 5–1,499, versus stated 6–1,500 | Make bounds and size universe explicit; do not claim identical selected pathways |
| BH is applied after selecting terms with CT overlap | `超幾何.txt:42`, 54, 84 | New workflow adjusts across the full size-eligible family, then filters displays |
| Hard-coded dimensions make scripts fragile | `超幾何.txt:29`, 79–81 use 7,857 columns and 100 rows | Use actual members and data dimensions |
| Existing saved score files do not match an exact negative-log transformation everywhere | Across 42,559 uniquely joinable genes, error >1e-6 occurs for 726 short and 13,807 middle values; max errors 0.0080 and 0.4897 | Record this observation; rounding, numerical cancellation or processing differences require tracing; do not overwrite the saved artifacts |
| Comments do not identify the fitted lambda reliably | Saved values imply median lambdas approximately 0.01 short and 0.47 middle; script comment says 0.99 | Store fitted parameters and full search; estimates from rounded saved values are not a rerun |
| Negative empty-index selection is unsafe | Import script lines 108/122 and comparison script lines 43–50 use `x[-which(...),]` | If `which` is empty, `-integer(0)` selects zero rows; use a logical complement |
| Plot labels can describe a different quantity or phase | Import script lines 129–143 logs values but labels the x-axis p-value; correction plotting edits short then plots middle | Use phase-specific named outputs and explicit `-log10` labels |
| Gene length and SNP-count objectives are mixed | Correction files cited above | Keep separate diagnostics; matching one does not imply matching the other |
| Bar-chart significance annotations are hard-coded | `條狀圖證明.txt:30–33`, 48–51; the middle/combined titles still say short networks | Derive annotations from saved tests; do not draw fixed lines as proof of significance |
| Group tests lack multiplicity adjustment and verified independence | `比較分數與無母分析程式.txt:98–188` | Define disjoint groups, a comparison family, effect sizes, and whether selection reused the same scores; no validation claim from these scripts alone |
| Crosstalk edge widths lose the edge-to-weight mapping | `crosstalk 20240916.R:47–61`, 89–100 apply `unique()` to weights, then index by edge order at 146–147 | Join edge attributes by canonical unordered gene-set pairs; use stored edge weights directly |
| Crosstalk figure depends on manual indices and an interactive GUI | Same file: `loc3=47` at 139, row selections at 150/161, `tkplot` at 163 | Compute metadata by pathway ID and save deterministic layouts and thresholds |
| Research metadata and runtime environment are implicit | Absolute paths, in-script `install.packages`, missing seed/session/version metadata | New workflow uses relative paths, base R, config snapshots, hashes, seeds, and session records |

The group-comparison and crosstalk scripts have been reviewed and archived, not completely rewritten. Network provenance, final group assignments, and the intended figure inputs need their own reproduction pass. Corrected correlation plots are diagnostic; they do not prove uniform null calibration or independent experimental validation.

## Decisions required for an exact paper reproduction

1. Define columns 2–6 of the external files, including Edward Kao's score formulas, directionality, and which score entered each method. The diagnostic runs check columns 5 and 6 separately and assume larger means stronger.
2. Confirm the scientific resolution of duplicated and ambiguous gene mappings. Averaging in the diagnostic run is explicit but provisional.
3. Reconcile the paper's FOSCO input with the recovered upstream p-value fields and confirm the SNP-window definition. The new imputation policy is specified and tested; the exceptional replacement and short score now have a numerical reconstruction hypothesis whose historical use remains unconfirmed.
4. Confirm the historical background, pathway boundaries, BH testing family, secondary selection rules, and final GSEA implementation. Current scripts and the described input files differ.
5. Clarify MaxMean's publication status. The author says its results were omitted; the final online article names it in methods and figure wording. The repository keeps it experimental without asserting that the article never mentioned it.
6. Compare an agreed full analysis against the final supplementary tables. This review does not conclude that any published biological finding is invalid; it identifies what must be resolved to reproduce the computations.

## Validation evidence

See `reports/VALIDATION.md` for completed numerical tests, real-data diagnostic runs, and limitations. See [METHODS.md](METHODS.md) for the implemented null models. In particular, gene-set sampling from the background is not the same as a self-contained phenotype-permutation test; see [Goeman and Bühlmann](https://pubmed.ncbi.nlm.nih.gov/17303618/).
