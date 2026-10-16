# Statistical definitions for the cleaned workflow

This file describes the new implementation, including deliberate departures from the historical scripts. Numerical results should not be presented as the original paper's results without reconciliation.

## Universe and pathway eligibility

Each phase uses the unique, resolved genes in its own external input as the measured universe. Its foreground is the intersection of that universe with the corresponding primary CTgene list. Missing CTgenes are reported, not assigned an invented external score. Each GO set is intersected with the same universe and deduplicated before its size is calculated.

The default size interval is **6 through 1,500 inclusive**, applied to measured set size. Historical scripts use `N.genes > 4 & N.genes < 1500`, i.e. 5 through 1,499 on the original annotation size. This difference is intentional and visible in `config.R`; the default follows the paper's stated bounds but does not recreate its historical selection procedure.

All eligible terms, including those with zero CTgenes, enter each method's BH adjustment, separately by phase and method. Candidate plots can subsequently display terms with CTgenes. This avoids defining the hypergeometric testing family by the observed overlap. It does not apply the paper's secondary count/ratio or top-fraction selection rules; those require a confirmed reproduction protocol.

## Hypergeometric test

Let N be universe size, K be CTgenes in that universe, n be measured genes in a pathway, and k their CTgene overlap. Calculate `phyper(k - 1, n, N - n, K, lower.tail = FALSE)` and adjust across the complete eligible family with BH. Direct upper-tail calculation avoids cancellation from subtracting a probability close to one. See the [R distribution documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html).

## Scores and random gene sets

The main configuration leaves the score choice unset. The supplied diagnostic and full configurations analyze source column 5 and column 6 separately as nonnegative, larger-is-stronger scores. This is an explicit analysis assumption, not a conclusion about how the collaborator constructed them. No primary CTgene score is added to an external score. The full configurations explicitly use the arithmetic mean when multiple rows resolve to one gene; the default policy otherwise rejects duplicates.

SUMSTAT sums scores for every measured gene in each pathway. Each null draw samples the same number of genes without replacement from the same phase's universe. Random draws may repeat across iterations. Null distributions are reused for pathways with the same measured size. This preserves each pathway's marginal null while reducing runtime.

GSEA sorts by decreasing score, breaking ties by gene ID, then calculates the maximum positive weighted hit-minus-miss running sum. Weights use exponent 1. A zero total hit weight falls back to equal weights. The gene-set null samples positions in this same ranked universe and recalculates their weights. This is a **positive-tail GSEA-style test with nominal p-values and BH adjustment**; it does not produce Broad GSEA's signed two-tail ES, NES, leading edge, or its empirical FDR statistic. Consult the [GSEA user guide](https://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/) when comparing implementations.

For SUMSTAT, GSEA, and MaxMean use `(b + 1) / (B + 1)`, where b counts null statistics greater than or equal to the observed statistic (with a small floating-point tolerance). Store the observed statistic, null exceedance count, seed, and B in every result row. This avoids zero p-values from finite random sampling; see [Phipson and Smyth](https://pubmed.ncbi.nlm.nih.gov/21044043/). At B = 199, the floor is 0.005; at B = 10,000 it is approximately 0.0001, before BH correction. A result at the floor does not establish a smaller p-value.

These implementations use **gene sampling against an external background**. Calling SUMSTAT or MaxMean self-contained merely because of the statistic's name would misdescribe this null. They do not use replicate-level phenotype permutations and do not preserve within-pathway gene correlations. See [Goeman and Bühlmann](https://pubmed.ncbi.nlm.nih.gov/17303618/). Likewise, an external dataset used to select pathways is not an untouched holdout for evaluating those same selected pathways.

## FOSCO

For a confirmed input p in (0,1] and positive SNP count M, calculate `-expm1(M^lambda * log1p(-p))`. This stably evaluates `1 - (1-p)^(M^lambda)`. Fit lambda separately for each phase by minimizing the absolute Spearman correlation of adjusted scores with SNP counts over 0.01, 0.02, ..., 1.00. Ties select the smallest grid value. Record the full search and before/after correlation.

The correction entry point rejects duplicate genes and zero/missing SNP counts. It does not silently bin-impute or change the target from SNP count to gene length. It does not automatically apply a correction to the supplied external scores.

The paper's fold-change-derived input `10^(-abs(log(FC)))` has `abs(log(FC))` as its negative-base-10-log score. Such a number is not automatically a calibrated statistical p-value. The generic `fold_change` score option performs this magnitude transformation without SNP correction; it discards direction. A uniform QQ appearance or a low correlation with size does not validate inferential calibration.

## Experimental MaxMean

Standardize the selected external scores across the entire universe. For each set calculate `max(mean(pmax(z,0)), mean(pmax(-z,0)))`, dividing both contributions by the complete set size. Use matched-size random gene sets and the same finite-sampling p-value rule. This does not add the extra within-pathway centering/scaling found in the January 27 script and is not presented as the full Efron–Tibshirani restandardization procedure.

MaxMean is excluded from default and 199-draw diagnostic runs. The full configurations explicitly include it to exercise the method on real data; it remains labelled experimental in the documentation and plots. Its experimental status follows the author's intended repository scope, while the paper's wording is recorded separately in the review.

## Matching Python and R implementations

`python/ctpathway.py` and `R/core.R` implement the same formulas. Python uses [SciPy's hypergeometric survival function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html), NumPy arrays, and an explicit BH implementation. R uses `phyper` and `p.adjust`. The cross-language reference test enumerates every three-gene set in a ten-gene universe, so its comparison of null statistics and p-values is independent of random generators.

Full simulations use R's native generator with `set.seed` and `sample.int`, and Python's `numpy.random.default_rng` with [sampling without replacement](https://numpy.org/doc/stable/reference/random/generated/numpy.random.Generator.choice.html). R uses hash sampling when the sampled set occupies at most half the universe to reduce allocation. A seed reproduces a run within its recorded implementation/environment; using the same seed in both languages does not produce identical random draws. Pathways sharing a size also share a simulated null within that method and phase, so their Monte Carlo errors are dependent.

The full-output verifier requires agreement of pathway membership, counts, observed statistics, and exact hypergeometric probabilities. It independently recalculates BH and checks every Monte Carlo p-value against its stored exceedance count. Differences between R and Python simulated p-values are reported with a broad six-standard-error diagnostic screen; this screen is not a formal proof of equivalence or a biological validation.
