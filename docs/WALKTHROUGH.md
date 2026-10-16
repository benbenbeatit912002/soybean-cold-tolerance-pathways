# How to read and run the R and Python versions

Both versions implement the same statistical definitions and read the same prepared input files. You can run the R analyses without Python, or the Python analyses without R. The cross-language validation test needs both interpreters. The original-data staging and provenance utilities are shared Python utilities; they are not required when using the prepared data already included in the package.

## Start with ten genes

Run either `Rscript examples/walkthrough.R` or `python examples/walkthrough.py`. The two short files explain the same example line by line. The example uses invented scores for ten genes, and a three-gene pathway with scores 5, 7, and 8. Its SUMSTAT is 20. Random three-gene sets show how unusual a sum of 20 is in that ten-gene background.

This example is for learning; its numbers are not soybean results.

## The real-data workflow

1. **Load the inputs.** The CTgene files define the primary candidates. The external files provide measured gene scores. The GO file defines which genes belong to each pathway. These three roles are separate.
2. **Use one gene universe per phase.** Keep each measured gene once when calculating membership. Restrict the CT list and every GO pathway to that same universe. Save the missing CTgenes and duplicate rows for inspection.
3. **Test overlap.** The hypergeometric test asks whether a pathway contains more CTgenes than expected from its size. It uses counts, not external score magnitudes.
4. **Test scores.** SUMSTAT adds all measured pathway scores. GSEA asks whether pathway genes concentrate toward the high-score end of the ranked universe. Experimental MaxMean compares positive and negative parts of globally standardized scores.
5. **Generate random sets.** For a pathway containing n measured genes, sample n distinct genes from the same universe and recalculate its statistic. Repeat 10,000 times in a full run. Genes can appear again in another draw. Pathways of the same size share a null distribution to save computation.
6. **Calculate p-values and adjust.** Count null statistics at least as large as the observed statistic; use `(count + 1)/(10000 + 1)`. Apply BH across the full size-eligible pathway family, separately for each method and phase.
7. **Write tables and figures.** The inputs, settings, duplicate policy, random seed, number of draws, raw p-values, and adjusted p-values remain inspectable.

## Find the corresponding code

| Purpose | R | Python |
|---|---|---|
| Main analysis | `scripts/run_analysis.R` | `scripts/run_analysis.py` |
| Statistical functions | `R/core.R` | `python/ctpathway.py` |
| Score correction | `scripts/correct_scores.R` | `scripts/correct_scores.py` |
| Plotting | `scripts/plot_results.R` | `scripts/plot_results.py` |
| Tiny teaching example | `examples/walkthrough.R` | `examples/walkthrough.py` |
| Full column-5 configuration | `config.score5-full.R` | `config.score5-full.json` |
| Full column-6 configuration | `config.score6-full.R` | `config.score6-full.json` |

Functions are kept small and named after their task. Main scripts use numbered comments so that a reader can follow the workflow without opening the archived exploratory scripts.

## Install and run

Use R 4.5.1 for the tested R environment; no additional R packages are needed. Use Python 3.11 or newer and install `requirements.txt` in your own environment:

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

Both correction scripts require `gene_id`, `raw_short`, `raw_middle`, and positive `snp` counts. The raw quantities must be in `(0,1]`; confirmed fold-change-derived pseudo-probabilities still must not be called calibrated p-values. Do not feed the supplied column-5/6 scores directly into FOSCO. Their score formulas and the intended raw input must first be established. Both implementations are tested on synthetic correction inputs; a real-data FOSCO rerun remains a separate provenance decision.

For the method assumptions and scope, read [METHODS.md](METHODS.md); for the original code findings, read [CODE_REVIEW.md](CODE_REVIEW.md).
