# Run the R workflow

The R implementation uses **base R only**. Python, RStudio, Excel, and additional R packages are not required. The tested version is R 4.5.1.

## 1. Open a terminal in the repository folder

Check your installation with `Rscript --version`. Then:

```sh
Rscript --vanilla run_pipeline.R --mode demo
```

If Rscript is not on PATH, use the installed executable in Windows PowerShell, adjusting the version:

```powershell
& "C:/Program Files/R/R-4.5.1/bin/Rscript.exe" --vanilla run_pipeline.R --mode demo
```

`--vanilla` starts a fresh R session without saved workspaces or personal startup files.

## 2. Check the demo

Find six method tables and six PNG figures in `results/r_demo/analysis/`. In the synthetic short-term SUMSTAT table, `TOY:1` has statistic **105**. These are teaching data, not soybean findings. The separate `Rscript --vanilla examples/walkthrough.R` explains a ten-gene example line by line.

## 3. Recalculate real data

Start with the shorter check, then run the complete reanalysis:

```sh
Rscript --vanilla run_pipeline.R --mode quick --use-reconstruction-hypotheses
Rscript --vanilla run_pipeline.R --mode full --use-reconstruction-hypotheses
```

For all implemented methods, including experimental MaxMean:

```sh
Rscript --vanilla run_pipeline.R --mode full --use-reconstruction-hypotheses --include-maxmean --output results/r_full_with_maxmean
```

Both real-data modes reconstruct scores from included expression, correspondence, length, and SNP marker files. Quick uses 199 draws per size/method; full uses 10,000. Hypergeometric probabilities need no simulations.

The hypotheses flag is a scientific choice. Read [reproducibility status](REPRODUCIBILITY.md). Omitting it saves the strict calculation and exits with status 2 when unresolved scores remain, before enrichment.

## 4. Read or repeat a run

`run_manifest.tsv` records mode, hypotheses, and completion. `execution_steps.tsv` records subprocess exit codes and durations. `analysis_config.R` records settings. Inspect the exception audit and input summary before interpretation.

Existing nonempty output folders are refused. For a repeat run, add `--output results/r_quick_repeat`.

The entry point locates the repository from its own path, so it works when called by full path from another directory. Defaults are anchored to the repository. Explicit relative output paths are relative to the caller's directory. Quote paths containing spaces.

## 5. Run R tests

```sh
Rscript --vanilla tests/test_core.R
Rscript --vanilla tests/test_integration.R
Rscript --vanilla tests/test_short_reference.R
Rscript --vanilla tests/test_historical_short.R
```

Python and R use different native random generators, so simulated p-values can differ while observed statistics agree. Optional cross-language checks require both interpreters; see the [Python guide](RUN_PYTHON.md#5-run-tests).

## Troubleshooting

| Message | Action |
|---|---|
| Rscript not found | Add R to PATH or invoke the quoted executable path |
| Output must be new or empty | Select another output directory |
| Upstream checksum mismatch | Restore packaged inputs; base R verifies their MD5 hashes |
| Exit status 2 | Read the unresolved score audit and reconstruction documentation |
| MaxMean warning | Expected when experimental MaxMean was requested |

Lower-level scripts remain in `scripts/`; run them from the repository root. `scripts/run_analysis.R` alone uses an older membership-only default. Use `run_pipeline.R` for the complete current workflow.
