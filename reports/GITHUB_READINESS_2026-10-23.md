# GitHub preparation and portable validation - 23 October 2026

## What changed

- Added complete independent Python and base-R entry points with demo, quick, and full modes.
- Included byte-preserved expression, gene-length, and SNP-marker inputs; reused the existing correspondence file. Relative manifest paths replace the original archive layout.
- Added SHA-256/MD5 input fingerprints, protected data bytes from Git line-ending conversion, and rejected mismatched sources.
- Added English README, separate Python/R guides, data provenance notes, reproducibility status, and packaging instructions.
- Retained original scripts and prior analyses as explicitly labelled archival/advanced material.
- Added fresh-directory packaging checks and configured Windows/Linux CI. Remote CI execution is not claimed.

## Complete recalculation and simulation check

A source ZIP with no saved results was extracted into a new directory containing spaces. Entry points were invoked from outside that repository. The included data supplied all calculation inputs.

The full checkpoint passed all 15 command groups, including both languages' tests, shared-null checks, demos, real-data quick runs, complete full runs, intentional strict-mode stops, and independent result verification. Its Python unit suite contained 24 tests. Two additional entry-point tests were added afterward and passed locally; the final package check reruns the expanded suite.

| Full run | Observed exit code | Elapsed seconds |
|---|---:|---:|
| Python | 0 | 637.67 |
| R | 0 | 1601.96 |

Both full runs include short and middle hypergeometric, SUMSTAT, positive-tail GSEA, and experimental MaxMean. Each simulation uses 10,000 draws per distinct measured pathway size. Across the two languages this is **45,240,000 random sets**, 12 simulation tables plus four exact hypergeometric tables, and 16 regenerated figures.

Each language recalculated 42,575 scores per phase. All short scores match the packaged archived reference at exact decimal precision. Both phases' exported scores agree between languages and with the previous 21 October recalculation. All 16 complete method tables match that earlier same-language run in every parsed field, including native Monte Carlo results with the same seed.

The independent verifier checks measured memberships, hypergeometric probabilities, SUMSTAT totals, draw counts, exceedance formulas, BH adjustment, cross-language observed statistics, and generated figure presence. Native R/Python Monte Carlo p-values are reported rather than forced equal.

Strict runs omitted the hypotheses and exited with the expected status 2, leaving one short and two middle nonfinite scores and no enrichment outputs. Default demo and no-overwrite behavior are covered by the added entry-point tests.

## Evidence and snapshots

- [Full execution checkpoint](portable_full_package_check_2026-10-23.json): tested ZIP and source-file hashes, commands, exit codes, durations, and retained logs.
- [Independent full-result verification](portable_full_results_2026-10-23.json).
- [Comparison with the previous recalculation](portable_previous_run_comparison_2026-10-23.json).
- [Latest package check](latest_portable_package_check.json): all 13 command groups passed, including the expanded 26-test Python suite, R tests, both real-data quick runs, expected strict-mode stops, and independent result verification.

Full-test ZIP SHA-256: `d9ceb3ef5bfcd51ad9c7856da2fb037def7ca0f95cc14d5244dedb6af819eba5`.

The full-test checkpoint precedes final documentation edits and two added tests. Analysis code, pipeline entry points, dependencies, and scientific data remain byte-identical to that checkpoint. Subsequent clean-package acceptance passed with the expanded test suite and both real-data quick runs. Unchanged full simulations were not repeated during documentation cleanup. The latest package report identifies the exact tested ZIP hash; this summary and the execution report were updated after that check. Because a report cannot include the hash of an archive containing that same report, consult the external/latest report for tested archive identity.

Plots inspected locally include Python full short SUMSTAT and R quick middle GSEA. The Windows R host prints inherited C.UTF-8 locale warnings and the expected experimental MaxMean warnings; the recorded commands pass despite those warnings.

## Scientific and publication limits

The score reconstruction hypotheses remain explicitly unconfirmed. The converted Yamasaki/Peixiu TXT is included; regenerating it from the Yamasaki Excel workbooks is not established. Published pathway reproduction remains unresolved; the packaging changes do not change that conclusion.

Software and data licensing remain unassigned. This report covers local testing with Windows, Python 3.11, and R 4.5.1. Remote Windows/Linux test outcomes are recorded separately by GitHub Actions.
