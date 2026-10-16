# Soybean cold tolerance pathways

Python and R workflows for short- and mid-term soybean cold-tolerance pathway analysis using hypergeometric testing, SUMSTAT, and GSEA.

Research companion for [Liu et al. (2025), *A multi-layered systems biology framework reveals dual-phased regulators and hormonal crosstalk underlying soybean cold tolerance*](https://doi.org/10.1007/s00299-025-03643-2).

This package organizes the short- and mid-term analysis, preserves the original scripts, and provides matching **Python and R workflows**. Start with the [step-by-step walkthrough](docs/WALKTHROUGH.md) to understand the code, then the [data dictionary](docs/DATA_DICTIONARY.md) and [code review](docs/CODE_REVIEW.md). This is a documented reanalysis; matching the published numerical results still requires confirmation of the historical score and analysis settings.

## Scope and inputs

- **Primary CTgenes:** `CTgenes 含分數/SHORT.txt` and `MIDDLE.txt`: 44 and 143 genes, respectively; 170 unique genes together. Their prior scores are preserved.
- **External data:** the two `WITHOUTINF/outer_data_2.0_*_address.txt` files specified by the author. Columns 5 and 6 are scores supplied by Edward Kao; their exact construction and historical usage remain to be documented.
- **GO annotation:** `GO_terms_4899.txt`, containing 4,896 GO terms despite the filename. The table contains GO IDs, descriptions, gene counts, and member genes. Gene scores are joined from separate inputs.
- **Methods:** hypergeometric, SUMSTAT, and a positive-tail weighted GSEA-style analysis. MaxMean is available as an experimental option.
- **Excluded:** long-term analyses, complete SoyNet/Cytoscape reconstruction, and experimental laboratory analyses. Related scripts in the requested code folder are reviewed and archived.

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
Rscript --vanilla scripts/run_analysis.R
Rscript --vanilla scripts/plot_results.R
```

On the author's Windows machine, the R executable is `C:/Program Files/R/R-4.5.1/bin/Rscript.exe`. In PowerShell, invoke that quoted path with `&` if `Rscript` is not on PATH.

The default run computes hypergeometric results for both phases, lists missing CTgenes and duplicate rows, and saves settings and environment information. It uses each external dataset's resolved gene IDs as its background. It does not require an external numerical score.

## Full simulations in both languages

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

These configurations assume larger supplied scores mean stronger evidence and explicitly average scores for duplicate gene IDs. These remain assumptions to confirm before interpreting the results scientifically. The minimum attainable nominal simulation p-value is `1/10001`, approximately 0.0001; BH correction is then applied across all eligible pathways.

The full runs and checks are recorded in [full simulation validation](reports/FULL_SIMULATION_VALIDATION.md). After all four runs finish, verify their tables and compare languages with:

```sh
python tests/check_cross_language.py
python scripts/verify_full_runs.py
```

Python and R agree on deterministic statistics and exact hypergeometric probabilities within numerical precision. Their native random-number generators produce different null samples, so simulated p-values can differ. The comparison checks identical enumerated null sets separately and reports the differences in the full runs.

For a quicker execution check, use `config.score5-diagnostic.json` / `.R` or their score-6 equivalents. These use only 199 draws and exclude MaxMean. The earlier R results in `score5_diagnostic` and `score6_diagnostic` are retained as historical diagnostics, not final significance estimates.

Original CT prior scores are not added to external scores. They are different measurements, and a combined score would require a defined formula. Short and middle are analyzed separately.

## Correction and plots

`scripts/correct_scores.py` and `scripts/correct_scores.R` apply numerically stable FOSCO separately to both phases. They require a table with `gene_id`, `raw_short`, `raw_middle`, and positive `snp` counts, with raw values in `(0,1]`. Each version records the lambda search and generates labelled diagnostics. Do not give these scripts column 5 or 6 directly: those supplied scores exceed one.

```sh
Rscript --vanilla scripts/correct_scores.R confirmed_raw_scores.tsv results/correction
python scripts/correct_scores.py confirmed_raw_scores.tsv results/python_correction
```

The correction entry point was tested with synthetic data. The historical correction was audited, but has not been reconstructed end to end from a confirmed raw-score definition. Removing a correlation does not establish that a transformed score is a calibrated p-value.

## Repository contents

| Location | Purpose |
|---|---|
| `R/core.R` and `python/ctpathway.py` | Matching statistical functions with explanatory comments |
| `config*.R` and `config*.json` | Explicit analysis and simulation settings |
| `examples/walkthrough.R` and `.py` | Matching ten-gene teaching examples |
| `docs/WALKTHROUGH.md` | Plain-English explanation of the complete workflow |
| `scripts/` | Input staging, audits, analysis, correction, and plotting |
| `data/input/` | Staged primary data, external scores, GO membership and ID correspondence |
| `legacy/` | Byte-identical copies of all 11 scripts from the requested code folder |
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

`python scripts/verify_outputs.py` checks the earlier membership/199-draw R outputs. `python scripts/verify_full_runs.py` checks the new full runs in both languages. `python scripts/package_project.py` creates a local review ZIP, including generated results and excluding temporary tests, virtual environments, and caches.

With that ZIP present, `python scripts/check_package.py` extracts it to a fresh folder and runs both languages' tests, examples, membership analyses, and plotting commands, then verifies the included full simulation outputs. Supply `--rscript "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"` if needed. This execution check does not repeat the expensive full simulations; it keeps its logs and writes `reports/latest_dual_language_package_check.json`.

## Publication status

Repository: [benbenbeatit912002/soybean-cold-tolerance-pathways](https://github.com/benbenbeatit912002/soybean-cold-tolerance-pathways). The repository includes code, prepared inputs, documentation, and validation reports. Generated `results/` files are excluded from Git and can be regenerated with the commands above; the author's local review ZIP also includes the completed results.

Before describing the analyses as an exact reproduction, resolve the decisions in the code review and compare a confirmed full run against the final supplementary tables. The published article mentions MaxMean; the author's intended repository treats its results as experimental. That distinction is documented rather than silently removed.

No software or data license has been assigned on the author's behalf. Select appropriate licenses and confirm redistribution terms for collaborator-provided data before public release. Cite the article using `CITATION.cff`.
