# Run the Python workflow

Python runs independently of R. Commands assume a terminal in the repository folder unless stated otherwise. The tested interpreter is Python 3.11; dependency versions are pinned in `requirements.txt`.

To include the historical Excel selection and offline crosstalk, add `--crosstalk reference` to a real-data command. See [crosstalk instructions](CROSSTALK.md) for complete commands and reuse of existing analysis results.

## 1. Install in a virtual environment

Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe run_pipeline.py --mode demo
```

Using the executable directly avoids PowerShell activation-policy changes. In subsequent Windows commands, replace `python` with `.\.venv\Scripts\python.exe` if needed.

macOS or Linux:

```sh
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python run_pipeline.py --mode demo
```

Dependencies are downloaded during installation. The analysis itself uses included files and needs no network connection.

## 2. Check the demo

Find six method tables and six PNG figures in `results/python_demo/analysis/`. In the synthetic short-term SUMSTAT table, `TOY:1` has statistic **105**. These are teaching data, not research findings. The separate `python examples/walkthrough.py` explains a ten-gene example line by line.

## 3. Recalculate real data

Run a shorter execution check, then the full analysis:

```sh
python run_pipeline.py --mode quick --use-reconstruction-hypotheses
python run_pipeline.py --mode full --use-reconstruction-hypotheses
```

For all implemented methods, including experimental MaxMean:

```sh
python run_pipeline.py --mode full --use-reconstruction-hypotheses --include-maxmean --output results/python_full_with_maxmean
```

Quick and full use the same complete inputs and freshly calculated scores. Quick uses 199 draws per size/method; full uses 10,000. Hypergeometric results are exact in either mode.

Read [reproducibility status](REPRODUCIBILITY.md) before enabling the hypotheses. To investigate the strict calculation, omit the flag and use a new output directory. Exit status 2 means unresolved values were saved and enrichment was not run.

## 4. Inspect or repeat

`run_manifest.json` records completion and subprocess exit codes. `analysis_config.json` captures every setting. Inspect `scores/exception_audit.tsv` and `analysis/input_summary.tsv` before interpreting results.

For a repeat run, add `--output results/python_quick_repeat`. Nonempty output folders are refused.

The entry point also works from another working directory when called by full path. Default paths are anchored to the repository. An explicit relative `--output` is relative to your terminal's current directory. Quote paths containing spaces.

## 5. Run tests

```sh
python -m unittest discover -s tests -p "test_*.py"
```

Tests requiring R are skipped if R is unavailable. If both languages are installed:

```sh
python tests/check_cross_language.py
python scripts/verify_portable_runs.py results/python_full results/r_full
```

The verifier needs completed real-data full runs in both languages with matching method choices. For quick runs, use their directories and add `--permutations 199`.

## Troubleshooting

| Message | Action |
|---|---|
| Python not found | Install Python or use its full executable path; some systems use `python3.11` |
| Missing dependency | Install requirements with the same interpreter used for analysis |
| Output must be new or empty | Choose another `--output` |
| Upstream checksum mismatch | Extract clean input files; do not resave byte-preserved sources in a spreadsheet editor |
| Exit status 2 | Inspect the unresolved score audit and reconstruction documentation |

Lower-level recalculation, analysis, and plotting scripts remain in `scripts/`; run them from the repository root. `scripts/run_analysis.py` alone uses an older membership-only default. The complete current entry point is `run_pipeline.py`.
