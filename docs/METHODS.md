# Statistical definitions for the cleaned workflow

This file describes the new implementation, including deliberate departures from the historical scripts. The historical short and middle scores are recalculated from upstream inputs with explicitly documented repair hypotheses. Saved scores are used only for independent verification; the former `--short-reference` copying option has been removed. Pathway results remain reanalyses unless reconciled with the paper.

## Universe and pathway eligibility

The portable main workflow uses the 42,567 unique genes in its freshly recalculated score table as the measured universe for each phase. The older WITHOUTINF configurations instead use the unique, resolved genes in each phase-specific external input. Its foreground is the intersection of that universe with the corresponding primary CTgene list. Missing CTgenes are reported, not assigned an invented external score. Each GO set is intersected with the same universe and deduplicated before its size is calculated.

The default size interval is **6 through 1,500 inclusive**, applied to measured set size. Historical scripts use `N.genes > 4 & N.genes < 1500`, i.e. 5 through 1,499 on the original annotation size. This difference is intentional and visible in `config.R`; the default follows the paper's stated bounds but does not recreate its historical selection procedure.

All eligible terms, including those with zero CTgenes, enter each method's BH adjustment, separately by phase and method. Candidate plots can subsequently display terms with CTgenes. This avoids defining the hypergeometric testing family by the observed overlap. It does not apply the paper's secondary count/ratio or top-fraction selection rules; those require a confirmed reproduction protocol.

## Hypergeometric test

Let N be universe size, K be CTgenes in that universe, n be measured genes in a pathway, and k their CTgene overlap. Calculate `phyper(k - 1, n, N - n, K, lower.tail = FALSE)` and adjust across the complete eligible family with BH. Direct upper-tail calculation avoids cancellation from subtracting a probability close to one. See the [R distribution documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html).

## Scores and random gene sets

The portable entry points select the named `score_short` and `score_middle` columns from fresh reconstruction, treating them as already logged, larger-is-stronger scores and averaging duplicate-gene scores. The older `config.python.json` and `config.R` defaults leave score choice unset and run membership tests only. The supplied diagnostic and full configurations analyze source column 5 and column 6 separately as nonnegative, larger-is-stronger scores. This is an explicit analysis assumption, not a conclusion about how the collaborator constructed them. No primary CTgene score is added to an external score. The full configurations explicitly use the arithmetic mean when multiple rows resolve to one gene; the default policy otherwise rejects duplicates.

SUMSTAT sums scores for every measured gene in each pathway. Each null draw samples the same number of genes without replacement from the same phase's universe. Random draws may repeat across iterations. Null distributions are reused for pathways with the same measured size. This preserves each pathway's marginal null while reducing runtime.

GSEA sorts by decreasing score, breaking ties by gene ID, then calculates the maximum positive weighted hit-minus-miss running sum. Weights use exponent 1. A zero total hit weight falls back to equal weights. The gene-set null samples positions in this same ranked universe and recalculates their weights. This is a **positive-tail GSEA-style test with nominal p-values and BH adjustment**; it does not produce Broad GSEA's signed two-tail ES, NES, leading edge, or its empirical FDR statistic. Consult the [GSEA user guide](https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/) when comparing implementations.

For SUMSTAT, GSEA, and MaxMean use `(b + 1) / (B + 1)`, where b counts null statistics greater than or equal to the observed statistic (with a small floating-point tolerance). Store the observed statistic, null exceedance count, seed, and B in every result row. This avoids zero p-values from finite random sampling; see [Phipson and Smyth](https://pubmed.ncbi.nlm.nih.gov/21044043/). At B = 199, the floor is 0.005; at B = 10,000 it is approximately 0.0001, before BH correction. A result at the floor does not establish a smaller p-value.

These implementations use **gene sampling against an external background**. Calling SUMSTAT or MaxMean self-contained merely because of the statistic's name would misdescribe this null. They do not use replicate-level phenotype permutations and do not preserve within-pathway gene correlations. See [Goeman and Bühlmann](https://pubmed.ncbi.nlm.nih.gov/17303618/). Likewise, an external dataset used to select pathways is not an untouched holdout for evaluating those same selected pathways.

## FOSCO

For a confirmed input p in (0,1] and positive SNP count M, calculate `-expm1(M^lambda * log1p(-p))`. This stably evaluates `1 - (1-p)^(M^lambda)`. Fit lambda separately for each phase by minimizing the absolute Spearman correlation of adjusted scores with SNP counts over 0.01, 0.02, ..., 1.00. Ties select the smallest grid value. Record the full search and before/after correlation.

Before fitting, the correction entry point replaces zero SNP counts with the arithmetic mean of positive input counts in the same 2,000-bp gene-length group. Groups are `[1,2001)`, `[2001,4001)`, and so on, including lengths above 100,000 bp. A group without positive counts uses the overall positive-input mean. Donors are fixed before any replacement, so imputed counts never become donors. Round to the nearest integer, with ties to even, matching R's `round()` and NumPy's `rint()`. Positive counts remain unchanged. This reanalysis policy was introduced on 20 October 2026; the original script included zeros in group means and lacked the global fallback.

SNP counts must be finite nonnegative integers. Missing/negative/fractional counts and an input with no positive donors are rejected. Positive integer `gene_length` is required when any count is zero. Imputation affects the counts used by both the FOSCO formula and its SNP-correlation objective; gene length only defines donor groups. Outputs preserve `snp_original`, record `snp_imputation_source`, and write each replacement, bin, donor count, and unrounded mean to `snp_imputation.tsv`.

Duplicate genes are rejected by default. The explicit `--duplicate-policy preserve` option retains all records for historical correction audits, including duplicate rows in the donor means and lambda fits. This is a record-level replay, not a choice to overweight repeated genes in downstream pathway sampling. The [real-data validation](../reports/SNP_IMPUTATION_VALIDATION_2026-10-20.md) uses this option to preserve 42,575 historical records representing 42,567 unique genes. The correction entry point does not automatically apply a correction to the separate supplied WITHOUTINF scores.

### Recalculate the historical short scores

The separate `reproduce_historical_short.py` and `.R` commands rebuild the mapped inputs and original imputation from the expression table, gene correspondence, lengths, and raw SNP blocks. They do not read saved correction inputs or expected scores. The original group mean includes zeros and has no fallback, leaving one zero-SNP record. The lambda search gives 0.01 for short scores. Ordinary subtraction is deliberately retained because the stable formula changes extremely small historical probabilities.

The calculators reproduce the observed export rules: 15-significant-digit probability serialization before `-log10`, zero logged scores replaced by `1e-16`, and the saved fixed/scientific decimal format. Those rules match 42,574 scores exactly. The default retains the one unresolved infinite score and returns a nonzero status.

An explicit `--zero-repair-hypothesis` tests a possible last-row repair, using positive whole-genome donors in the same bin, the fitted middle-phase lambda, and seven-significant-digit probability rounding. It calculates the final saved value without loading it. Both languages then match all 42,575 exported short scores, but this hypothesis is not confirmed as the original procedure and is not a recommended new statistical method. See the [derivation and validation report](../reports/SHORT_SCORE_RECALCULATION_2026-10-21.md). The independent verifier checks exact decimals and input identity separately.

The old score-copying option was removed. Archived preservation reports document copying only and cannot serve as recalculation evidence.

The paper's fold-change-derived input `10^(-abs(log(FC)))` has `abs(log(FC))` as its negative-base-10-log score. Such a number is not automatically a calibrated statistical p-value. The generic `fold_change` score option performs this magnitude transformation without SNP correction; it discards direction. A uniform QQ appearance or a low correlation with size does not validate inferential calibration.

## Experimental MaxMean

Standardize the selected external scores across the entire universe. For each set calculate `max(mean(pmax(z,0)), mean(pmax(-z,0)))`, dividing both contributions by the complete set size. Use matched-size random gene sets and the same finite-sampling p-value rule. This does not add the extra within-pathway centering/scaling found in the January 27 script and is not presented as the full Efron–Tibshirani restandardization procedure.

MaxMean is excluded from default and 199-draw diagnostic runs. The full configurations explicitly include it to exercise the method on real data; it remains labelled experimental in the documentation and plots. Its experimental status defines the repository scope; the paper's wording is recorded separately in the review.

## Matching Python and R implementations

`python/ctpathway.py` and `R/core.R` implement the same formulas. Python uses [SciPy's hypergeometric survival function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html), NumPy arrays, and an explicit BH implementation. R uses `phyper` and `p.adjust`. The cross-language reference test enumerates every three-gene set in a ten-gene universe, so its comparison of null statistics and p-values is independent of random generators.

Full simulations use R's native generator with `set.seed` and `sample.int`, and Python's `numpy.random.default_rng` with [sampling without replacement](https://numpy.org/doc/stable/reference/random/generated/numpy.random.Generator.choice.html). R uses hash sampling when the sampled set occupies at most half the universe to reduce allocation. A seed reproduces a run within its recorded implementation/environment; using the same seed in both languages does not produce identical random draws. Pathways sharing a size also share a simulated null within that method and phase, so their Monte Carlo errors are dependent.

The full-output verifier requires agreement of pathway membership, counts, observed statistics, and exact hypergeometric probabilities. It independently recalculates BH and checks every Monte Carlo p-value against its stored exceedance count. Differences between R and Python simulated p-values are reported with a broad six-standard-error diagnostic screen; this screen is not a formal proof of equivalence or a biological validation.

## Recalculated-score run, 21 October 2026

The current recalculated-score workflow and its explicit settings are documented in
[the full rerun report](../reports/RECALCULATED_ANALYSIS_2026-10-21.md). It rebuilds
42,575 records from four upstream inputs, then tests 42,567 unique genes after
averaging duplicate scores. Both phases are run independently in Python and R.
The repair hypotheses, all-family BH analysis, historical SUMSTAT sensitivity
run, and intentionally invalid hypergeometric diagnostic profiles must remain
distinguishable. Exact saved-score agreement does not establish paper reproduction.
