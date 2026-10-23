# How to read and run the R and Python versions

For current installation and complete execution commands, use the [Python guide](RUN_PYTHON.md) or [R guide](RUN_R.md). The recommended entry points are `run_pipeline.py` and `run_pipeline.R`, using included data. Older branch-specific commands below explain previous analyses; see [reproducibility status](REPRODUCIBILITY.md).

Both versions implement the same statistical definitions and read the same prepared input files. You can run the R analyses without Python, or the Python analyses without R. The cross-language validation test needs both interpreters. The original-data staging and provenance utilities are shared Python utilities; they are not required when using the prepared data already included in the package.

## Start with ten genes

Run either `Rscript examples/walkthrough.R` or `python examples/walkthrough.py`. The two short files explain the same example line by line. The example uses invented scores for ten genes, and a three-gene pathway with scores 5, 7, and 8. Its SUMSTAT is 20. Random three-gene sets show how unusual a sum of 20 is in that ten-gene background.

This example is for learning; its numbers are not soybean results.

## The real-data workflow

1. **Load the inputs.** The CTgene files define the primary candidates. The packaged upstream files are used to recalculate measured gene scores. The GO file defines which genes belong to each pathway. These three roles are separate.
2. **Use one gene universe per phase.** Keep each measured gene once when calculating membership. Restrict the CT list and every GO pathway to that same universe. Save the missing CTgenes and duplicate rows for inspection.
3. **Test overlap.** The hypergeometric test asks whether a pathway contains more CTgenes than expected from its size. It uses counts, not external score magnitudes.
4. **Test scores.** SUMSTAT adds all measured pathway scores. GSEA asks whether pathway genes concentrate toward the high-score end of the ranked universe. Experimental MaxMean compares positive and negative parts of globally standardized scores.
5. **Generate random sets.** For a pathway containing n measured genes, sample n distinct genes from the same universe and recalculate its statistic. Repeat 10,000 times in a full run. Genes can appear again in another draw. Pathways of the same size share a null distribution to save computation.
6. **Calculate p-values and adjust.** Count null statistics at least as large as the observed statistic; use `(count + 1)/(10000 + 1)`. Apply BH across the full size-eligible pathway family, separately for each method and phase.
7. **Write tables and figures.** The inputs, settings, duplicate policy, random seed, number of draws, raw p-values, and adjusted p-values remain inspectable.

## Find the corresponding code

| Purpose | R | Python |
|---|---|---|
| Lower-level pathway analysis | `scripts/run_analysis.R` | `scripts/run_analysis.py` |
| Complete workflow | `run_pipeline.R` | `run_pipeline.py` |
| Historical score reconstruction | `scripts/recalculate_phase_scores.R` | `scripts/recalculate_phase_scores.py` |
| Statistical functions | `R/core.R` | `python/ctpathway.py` |
| Score correction | `scripts/correct_scores.R` | `scripts/correct_scores.py` |
| Plotting | `scripts/plot_results.R` | `scripts/plot_results.py` |
| Tiny teaching example | `examples/walkthrough.R` | `examples/walkthrough.py` |
| Full column-5 configuration | `config.score5-full.R` | `config.score5-full.json` |
| Full column-6 configuration | `config.score6-full.R` | `config.score6-full.json` |

Functions are kept small and named after their task. Main scripts use numbered comments so that a reader can follow the workflow without opening the archived exploratory scripts.

## Earlier WITHOUTINF branch (advanced)

For the current main workflow, follow the language guides above. The older commands in this section analyze supplied WITHOUTINF columns and do not recalculate the historical scores. The tested environments are R 4.5.1 and Python 3.11:

```sh
python -m pip install -r requirements.txt
```

Run commands from the repository root. For the full simulations in R:

```sh
Rscript --vanilla scripts/run_analysis.R config.score5-full.R
Rscript --vanilla scripts/run_analysis.R config.score6-full.R
Rscript --vanilla scripts/plot_results.R results/r_score5_full
Rscript --vanilla scripts/plot_results.R results/r_score6_full
```

For the matching Python workflow:

```sh
python scripts/run_analysis.py config.score5-full.json
python scripts/run_analysis.py config.score6-full.json
python scripts/plot_results.py results/python_score5_full
python scripts/plot_results.py results/python_score6_full
```

Each configuration covers both short and middle phases. It selects 10,000 draws and the three score methods, including experimental MaxMean. Hypergeometric results are calculated too. Do not run multiple copies writing to the same output directory. Use distinct `output_dir` settings for alternative analyses.

The full configurations retain the explicit diagnostic assumptions: score columns are analyzed separately as supplied, larger values mean stronger evidence, and repeated scores for a gene are averaged. These choices have not been established as the historical paper settings. A completed full simulation does not settle score provenance or establish paper reproduction.

## Read a result table

| Column | Meaning |
|---|---|
| `go_id`, `description` | Pathway identifier and name |
| `annotation_size` | Unique genes in the original GO annotation |
| `tested_size` | Pathway genes actually present in the phase's external universe |
| `ct_count`, `ct_ratio` | Measured CTgenes in the pathway and their fraction |
| `statistic` | Observed SUMSTAT, positive-tail GSEA-style ES, or experimental MaxMean |
| `null_exceedances` | Number of random sets at least as extreme, including ties |
| `permutations`, `seed` | Draw count and starting random seed |
| `p_value` | `(null_exceedances + 1)/(permutations + 1)` |
| `p_bh` | BH-adjusted p-value across all eligible terms |

For the exact hypergeometric test, the last two columns are instead named `hypergeometric_p` and `hypergeometric_bh`; it does not need simulation.

Some plotted bars can have exactly the same height. Finite simulations impose a minimum p-value, and BH adjustment can also produce ties. Equal bars therefore do not mean that the pathways have identical scores or biological effects. Consult the table's observed statistic, exceedance count, and unadjusted p-value alongside the adjusted value.

## Why two languages may give slightly different p-values

The observed statistics, gene counts, and exact hypergeometric probabilities should agree to numerical precision. R and NumPy use different random generators and sampling implementations. The same seed therefore does not mean the same random gene sets. Their Monte Carlo p-values can differ within sampling uncertainty, especially near a decision threshold.

`tests/check_cross_language.py` compares all three statistics using exactly the same 120 possible three-gene subsets of a ten-gene universe, and checks fixed and fitted FOSCO corrections. That test checks the formulas independently of RNG differences. Full-run comparison separately reports native Monte Carlo differences; it does not force results to match by copying one language's p-values.

## Correction is a separate, explicit step

Both correction scripts require `gene_id`, `raw_short`, `raw_middle`, and nonnegative integer `snp` counts. Include `gene_length` in base pairs if any SNP count is zero. For each zero, use the mean of nonzero SNP counts in the same 2,000-bp length group; if none exist, use the overall nonzero mean. Round the result to an integer, with ties to even. Only original positive counts are donors, so row order does not change the replacements.

For example, donor counts of 2 and 5 in one group give a replacement of 4. If another group has no nonzero counts, its replacement uses the positive counts across the entire input. Neither zero counts nor newly imputed counts enter these means. The original script included zeros and did not provide this global fallback.

```sh
python scripts/correct_scores.py raw_scores.tsv results/python_correction
Rscript --vanilla scripts/correct_scores.R raw_scores.tsv results/r_correction
```

Inspect `snp_imputation.tsv` for each replacement and its donor mean, `correction_summary.tsv` for the selected lambda, and `corrected_scores.tsv` for both original and imputed SNP counts plus phase-specific scores. Duplicate genes are rejected by default; `--duplicate-policy preserve` is available for an explicit historical record-level audit.

The raw quantities must be in `(0,1]`; fold-change-derived pseudo-probabilities still must not be called calibrated p-values. Do not feed the supplied column-5/6 scores directly into FOSCO. Both implementations passed synthetic tests and a [real-data rerun](../reports/SNP_IMPUTATION_VALIDATION_2026-10-20.md) from the saved historical input. That input already contains earlier SNP imputations. The WITHOUTINF configurations remain separate; the new `config.recalculated-full` configurations explicitly connect fresh historical-branch scores to pathway analysis.

## Recalculate and compare the original short scores

Use `scripts/reproduce_historical_short.py` or `.R` with `data/upstream` and an output directory (the original project layout is also supported). These readable entry points call `python/historical_short.py` and `R/historical_short.R`, respectively. Each implementation independently rebuilds the data from four upstream files; saved scores are never calculation inputs.

Follow the strict and optional hypothesis commands in the [advanced workflow reference](ADVANCED_WORKFLOWS.md#recalculate-the-historical-short-scores). Inspect `historical_snp_imputation.tsv`, `short_lambda_grid.tsv`, and `recalculated_short.tsv`. The latter preserves raw source decimals, raw marker counts, historical imputed counts, final counts, corrected probabilities, unrounded logarithms, and exported scores.

The strict run leaves one unresolved score and returns status 2 after writing its outputs. The optional hypothesis additionally writes `zero_repair_hypothesis.tsv` with actual donor IDs, their counts, the mean, fitted middle lambda, and probability rounding. The independent `verify_short_recalculation.py` command compares both versions in both languages against original files. It does not change any computed result.

The [report](../reports/SHORT_SCORE_RECALCULATION_2026-10-21.md) distinguishes the 42,574 values directly reconstructed with the export rules from the last value requiring an unconfirmed repair. The optional hypothesis matches all saved decimal values, but numerical equality does not prove historical provenance. The former `--short-reference` option that copied scores is no longer accepted.

For the method assumptions and scope, read [METHODS.md](METHODS.md); for the original code findings, read [CODE_REVIEW.md](CODE_REVIEW.md).

## Run the recalculated-score analysis

Follow [the two-phase commands in the README](../README.md#run-the-soybean-analysis).
The sequence is `recalculate_phase_scores` -> `run_analysis` -> `plot_results`.
Both Python and R implementations perform the same stages independently.

First inspect the score file's `status_short` and `status_middle` columns and
`exception_audit.tsv`. Then inspect `input_summary.tsv` in the analysis folder:
it separates input records, unique genes, measured CTgenes, and tested pathways.
The numerical score table is `short_analysis_scores.tsv` or
`middle_analysis_scores.tsv`. Each method then has one result table per phase.

In a simulation table, `statistic` is the observed pathway value,
`null_exceedances` counts how many of 10,000 random sets meet or exceed it, and
`p_value` is `(null_exceedances+1)/10001`. `p_bh` adjusts all eligible pathways
together. A zero exceedance count therefore has a finite p-value near 0.0001.
Hypergeometric probabilities are calculated directly and need no random draws.

The historical diagnostic command is separate because it preserves the
archived record-level sampling and `b/10000` convention. Its zero values and
intentionally invalid hypergeometric variants must not be substituted into
the main analysis. The [comparison report](../reports/RECALCULATED_ANALYSIS_2026-10-21.md)
explains which manuscript values remain inconsistent.
