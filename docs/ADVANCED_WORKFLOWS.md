> Historical reference, retained from the pre-portability README. New users should start with [the current README](../README.md). Commands below that name the original project directory are maintainer audits, not prerequisites for the portable workflow. The current packaging commands are in [RELEASING.md](RELEASING.md).

# Soybean cold tolerance pathways

Python and R workflows for short- and mid-term soybean cold-tolerance pathway analysis using hypergeometric testing, SUMSTAT, and GSEA.

Research companion for [Liu et al. (2025), *A multi-layered systems biology framework reveals dual-phased regulators and hormonal crosstalk underlying soybean cold tolerance*](https://doi.org/10.1007/s00299-025-03643-2).

This package organizes the short- and mid-term analysis, preserves the original scripts, and provides matching **Python and R workflows**. Start with the [step-by-step walkthrough](WALKTHROUGH.md) to understand the code, then the [data dictionary](DATA_DICTIONARY.md) and [code review](CODE_REVIEW.md). This is a documented reanalysis; matching the published numerical results still requires confirmation of the historical score and analysis settings.

**Recalculated-score analysis, 21 October 2026:** the new workflow below calculates both phases from upstream data and passes the resulting named score columns directly into the pathway analysis. Read the [full rerun and manuscript comparison](../reports/RECALCULATED_ANALYSIS_2026-10-21.md) for execution results, the two explicitly unconfirmed repair hypotheses, and differences from the manuscript. This workflow uses the historical expression-p-value branch, separate from the earlier WITHOUTINF column analyses.

**Short-score recalculation, 21 October 2026:** Python and R now rebuild the inputs from four upstream files and calculate the short scores. The archived procedure plus reconstructed export precision matches 42,574 saved values exactly. An optional, explicitly unconfirmed zero-SNP repair hypothesis also calculates the last value, producing 42,575 exact matches. The repair's historical provenance remains unresolved. See the [recalculation report](../reports/SHORT_SCORE_RECALCULATION_2026-10-21.md). Saved scores are used only for independent verification; the previous score-copying option has been removed.

**Paper concordance, 18 October 2026:** the [audit](../reports/PAPER_CONCORDANCE_2026-10-18.md) identifies differences in score preparation, gene universes, pathway eligibility, and final selection. Execution tests passed, but paper reproduction is not verified. In particular, column 5 is fold change; the existing raw-column-5 simulations do not apply the reported score preparation.

The [historical data-flow audit](../reports/HISTORICAL_DATA_FLOW_2026-10-18.md) reconstructs the original mapping and both saved FOSCO transformations, distinguishes record counts from unique genes, and links a data-dependency inventory of all 11 original scripts.

**SNP correction update, 20 October 2026:** both languages now replace zero SNP counts with the same length group's nonzero mean, falling back to the overall nonzero mean. A [real-data rerun and comparison](../reports/SNP_IMPUTATION_VALIDATION_2026-10-20.md) checked all 42,575 saved records in both phases. This new policy differs from the original script, whose means included zeros and had no global fallback.

## Scope and inputs

- **Primary CTgenes:** `CTgenes 含分數/SHORT.txt` and `MIDDLE.txt`: 44 and 143 genes, respectively; 170 unique genes together. Their prior scores are preserved.
- **External data:** the two `WITHOUTINF/outer_data_2.0_*_address.txt` files. Historical source matching identifies column 5 as fold change. Column 6 was supplied by Edward Kao and is historically labelled `Pvalue`, but exceeds one and has an unresolved transformation. See the [data dictionary](DATA_DICTIONARY.md).
- **GO annotation:** `GO_terms_4899.txt`, containing 4,896 GO terms despite the filename. The table contains GO IDs, descriptions, gene counts, and member genes. Gene scores are joined from separate inputs.
- **Methods:** hypergeometric, SUMSTAT, and a positive-tail weighted GSEA-style analysis. MaxMean is available as an experimental option.
- **Excluded:** long-term analyses, complete SoyNet/Cytoscape reconstruction, and experimental laboratory analyses. Related scripts in the original code archive are reviewed and archived.

## Quick start

Clone the repository and enter its folder:

```sh
git clone https://github.com/benbenbeatit912002/soybean-cold-tolerance-pathways.git
cd soybean-cold-tolerance-pathways
```

Run commands from this folder. Choose either language for analysis, correction, and plotting. The staged inputs are already included. R uses base R only. Python uses NumPy, SciPy, and Matplotlib; the tested versions are pinned in `requirements.txt`. Local validation used R 4.5.1 and Python 3.11.

**Python**

```sh
python -m pip install -r requirements.txt
python examples/walkthrough.py
python -m unittest discover -s tests -p "test_*.py"
python scripts/run_analysis.py
python scripts/plot_results.py
```

**R**

```sh
Rscript --vanilla examples/walkthrough.R
Rscript --vanilla tests/test_core.R
Rscript --vanilla tests/test_integration.R
Rscript --vanilla tests/test_short_reference.R
Rscript --vanilla tests/test_historical_short.R
Rscript --vanilla scripts/run_analysis.R
Rscript --vanilla scripts/plot_results.R
```

For example, a Windows R installation may place the executable at `C:/Program Files/R/R-4.5.1/bin/Rscript.exe`. In PowerShell, invoke that quoted path with `&` if `Rscript` is not on PATH.

The default run computes hypergeometric results for both phases, lists missing CTgenes and duplicate rows, and saves settings and environment information. It uses each external dataset's resolved gene IDs as its background. It does not require an external numerical score.

## Recalculate both phases and run full simulations

The packaged data/upstream manifest supplies four inputs: expression measurements, gene-ID correspondence, gene lengths, and SNP marker blocks. They are fingerprinted by the recalculation command. Saved gene scores are read only by verification commands. Run the following commands from the repository root.

```sh
python scripts/recalculate_phase_scores.py data/upstream results/recalculated_2026-10-21/python_scores --use-reconstruction-hypotheses
python scripts/run_analysis.py config.recalculated-full.json
python scripts/plot_results.py results/recalculated_2026-10-21/python_full
```

```sh
Rscript --vanilla scripts/recalculate_phase_scores.R data/upstream results/recalculated_2026-10-21/r_scores --use-reconstruction-hypotheses
Rscript --vanilla scripts/run_analysis.R config.recalculated-full.R
Rscript --vanilla scripts/plot_results.R results/recalculated_2026-10-21/r_full
```

The explicit hypotheses flag enables the documented zero-SNP repair and a middle-phase infinite-score replacement calculated from the largest finite score. Numerical agreement with the saved scores does not prove these were the historical manual operations. Without the flag, unresolved nonfinite scores are reported with exit status 2; they cannot enter enrichment analysis.

Each language tests both phases with exact hypergeometric probabilities and 10,000 random sets **per distinct measured pathway size and method** for SUMSTAT, positive-tail weighted GSEA, and experimental MaxMean. The analysis uses 42,567 unique measured genes, averages the eight extra duplicate records, tests measured pathway sizes 6 through 1,500, and applies BH across the entire eligible family. It does not automatically reproduce the older manuscript's testing family or GSEA implementation.

An additional historical SUMSTAT replay retains all 42,575 records, selects annotation sizes 5 through 1,499 with a CTgene overlap, and uses the archived strict-exceedance `b/B` convention. It also calculates explicitly labeled diagnostic hypergeometric variants, including the historical matrix-counting mistakes. Those variants are for explaining old results, not recommended statistical practice.

```sh
python scripts/run_historical_diagnostics.py results/recalculated_2026-10-21/python_scores/recalculated_scores.tsv results/recalculated_2026-10-21/python_historical
Rscript --vanilla scripts/run_historical_diagnostics.R results/recalculated_2026-10-21/r_scores/recalculated_scores.tsv results/recalculated_2026-10-21/r_historical
```

The [report](../reports/RECALCULATED_ANALYSIS_2026-10-21.md) gives the local manuscript extraction and independent verification commands. The verifier checks both phases against archived scores, every simulation's draw count and p-value formula, independent BH adjustment, cross-language statistics, and manuscript tables. Full outputs are generated under `results/recalculated_2026-10-21/` and are not tracked by Git.

## Earlier full simulations using WITHOUTINF columns

Each full configuration runs SUMSTAT, positive-tail GSEA, and experimental MaxMean with **10,000 random gene sets per distinct pathway size**, for both phases. Edward Kao's two score columns are analyzed separately. Hypergeometric tests are also saved; they are calculated directly and require no simulations.

```sh
python scripts/run_analysis.py config.score5-full.json
python scripts/run_analysis.py config.score6-full.json
python scripts/plot_results.py results/python_score5_full
python scripts/plot_results.py results/python_score6_full
```

```sh
Rscript --vanilla scripts/run_analysis.R config.score5-full.R
Rscript --vanilla scripts/run_analysis.R config.score6-full.R
Rscript --vanilla scripts/plot_results.R results/r_score5_full
Rscript --vanilla scripts/plot_results.R results/r_score6_full
```

These configurations assume larger supplied scores mean stronger evidence and explicitly average scores for duplicate gene IDs. The concordance audit established that this treats column 5 as raw fold change, without reciprocal-symmetric transformation or real-data FOSCO; those runs are execution records, not paper-reproduction results. Column 6's construction and the scientific duplicate policy remain unresolved. The minimum attainable nominal simulation p-value is `1/10001`, approximately 0.0001; BH correction is then applied across all eligible pathways.

The full runs and checks are recorded in [full simulation validation](../reports/FULL_SIMULATION_VALIDATION.md). After all four runs finish, verify their tables and compare languages with:

```sh
python tests/check_cross_language.py
python scripts/verify_full_runs.py
```

Python and R agree on deterministic statistics and exact hypergeometric probabilities within numerical precision. Their native random-number generators produce different null samples, so simulated p-values can differ. The comparison checks identical enumerated null sets separately and reports the differences in the full runs.

For a quicker execution check, use `config.score5-diagnostic.json` / `.R` or their score-6 equivalents. These use only 199 draws and exclude MaxMean. The earlier R results in `score5_diagnostic` and `score6_diagnostic` are retained as historical diagnostics, not final significance estimates.

Original CT prior scores are not added to external scores. They are different measurements, and a combined score would require a defined formula. Short and middle are analyzed separately.

## Correction and plots

### Recalculate the historical short scores

The historical calculators read the original expression table, correspondence table, gene lengths, and SNP marker blocks. They accept the packaged data/upstream manifest or the original project directory. They reconstruct the old imputation, search all 100 lambda values, and calculate and export the short scores. Both preserve source row identity and duplicate records.

```sh
python scripts/reproduce_historical_short.py data/upstream results/short_derivation_2026-10-21/python
Rscript --vanilla scripts/reproduce_historical_short.R data/upstream results/short_derivation_2026-10-21/r
```

The archived procedure leaves `Glyma.14G129500` at SNP zero. These commands write the outputs and return status 2 to report its unresolved infinite score; the other 42,574 exported scores match the original decimal values exactly. The observed output rules include intermediate probability rounding and replacement of zero logged scores by `1e-16`.

The [report](../reports/SHORT_SCORE_RECALCULATION_2026-10-21.md) explains an optional repair hypothesis: whole-genome positive SNP donors in the same length bin have mean 11.5, which rounds to 12; applying the fitted middle lambda to the remaining zero row and rounding its probability to seven significant digits produces the last saved short score. This is numerically consistent with the files but the historical operation is unconfirmed.

```sh
python scripts/reproduce_historical_short.py data/upstream results/short_derivation_2026-10-21/python_hypothesis --zero-repair-hypothesis
Rscript --vanilla scripts/reproduce_historical_short.R data/upstream results/short_derivation_2026-10-21/r_hypothesis --zero-repair-hypothesis
python scripts/verify_short_recalculation.py "D:/Bioinformatics project" results/short_derivation_2026-10-21
```

The archival verify_short_recalculation.py command above requires the original project directory and reads saved scores only after calculation. For portable verification, use scripts/verify_portable_runs.py as described in the current execution guides. It compares all 42,575 values using exact decimal equality, with no score tolerance. The hypothesis run matches every score in both languages, but its report retains the unresolved provenance flag. These commands do not establish agreement between the archived p-value input branch and the paper's fold-change description.

### Recompute FOSCO scores

`scripts/correct_scores.py` and `scripts/correct_scores.R` apply numerically stable FOSCO separately to both phases. They require a table with `gene_id`, `raw_short`, `raw_middle`, and nonnegative integer `snp` counts, with raw values in `(0,1]`. Include positive integer `gene_length` in base pairs whenever any SNP count is zero.

Zero SNP counts use the nonzero mean within 2,000-bp length groups (`[1,2001)`, `[2001,4001)`, and so on), or the overall nonzero mean if the group has no donors. Means use the original positive input counts and are rounded to the nearest integer, with ties to even. Both versions preserve original counts, write a row-level imputation audit, record the lambda search, and generate labelled diagnostics. Do not give these scripts column 5 or 6 directly: those supplied scores exceed one.

```sh
Rscript --vanilla scripts/correct_scores.R confirmed_raw_scores.tsv results/correction
python scripts/correct_scores.py confirmed_raw_scores.tsv results/python_correction
```

The correction entry points passed synthetic tests and a complete rerun from the saved historical `gene_p_len.txt`. Its remaining zero was replaced by 11, and recomputed Python/R scores agree within `5e-14`. That saved input already contains earlier imputations; this is not a reconstruction from raw SNP markers. The [20 October report](../reports/SNP_IMPUTATION_VALIDATION_2026-10-20.md) records that reanalysis; use the historical recalculation commands above to investigate the saved short scores. The older full configurations use the separate WITHOUTINF branch. The new `config.recalculated-full` configurations explicitly consume freshly rebuilt historical scores, including their documented repair hypotheses. Removing a correlation does not establish that a transformed score is a calibrated p-value.

## Repository contents

| Location | Purpose |
|---|---|
| `R/core.R` and `python/ctpathway.py` | Matching statistical functions with explanatory comments |
| `R/historical_short.R` and `python/historical_short.py` | Independent reconstruction and calculation of historical short scores |
| `config*.R` and `config*.json` | Explicit analysis and simulation settings |
| `examples/walkthrough.R` and `.py` | Matching ten-gene teaching examples |
| `docs/WALKTHROUGH.md` | Plain-English explanation of the complete workflow |
| `scripts/` | Input staging, audits, analysis, correction, and plotting |
| `data/input/` | Staged primary data, external scores, GO membership and ID correspondence |
| `data/reference/` | Expected historical inputs and scores for independent verification, with source hashes and record identities |
| `legacy/` | Byte-identical copies of all 11 scripts from the original code archive |
| `docs/CODE_REVIEW.md` | Findings, affected lines, consequences, and remediation |
| `docs/METHODS.md` | Statistical definitions and differences from historical scripts |
| `reports/` | Hash manifest, input audit, import decisions, parse checks, and validation summary |
| `results/` | Generated outputs, excluded from Git by default |
| `tests/` | Numerical checks and synthetic integration checks |

Regenerate the data package from the original folder:

```sh
python scripts/audit_inputs.py "D:/Bioinformatics project"
python scripts/prepare_inputs.py "D:/Bioinformatics project"
python scripts/audit_lineage.py "D:/Bioinformatics project"
Rscript --vanilla scripts/check_legacy.R
```

The originals are unchanged. The manifest records SHA-256 hashes and every staging transformation. Import-decision tables record every external row, including exclusions. Old-format IDs are mapped only when the supplied correspondence table has exactly one target.

`python scripts/verify_outputs.py` checks the earlier membership/199-draw R outputs. `python scripts/verify_full_runs.py` checks the new full runs in both languages. Current source and optional results packaging commands are documented in [RELEASING.md](RELEASING.md).

The current `check_package.py` tests a source-only ZIP from a clean extraction. Its `--full` option reruns complete simulations. See [RELEASING.md](RELEASING.md) for current commands; the older package-check report remains an archival execution record.

## Publication status

Repository: [benbenbeatit912002/soybean-cold-tolerance-pathways](https://github.com/benbenbeatit912002/soybean-cold-tolerance-pathways). The repository includes code, prepared inputs, documentation, and validation reports. Generated `results/` files are excluded from Git and can be regenerated with the commands above. An optional results bundle can be built using RELEASING.md.

Before describing the analyses as an exact reproduction, resolve the decisions in the code review and compare a confirmed full run against the final supplementary tables. The published article mentions MaxMean; this repository treats its results as experimental.

Software and data licenses remain unassigned. Select appropriate licenses and confirm redistribution terms for collaborator-provided data before public release. Cite the article using `CITATION.cff`.
