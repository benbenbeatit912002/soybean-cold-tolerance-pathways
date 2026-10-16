# Python and R full simulation validation — 15 October 2026

**All four full configurations completed successfully.** SUMSTAT, positive-tail GSEA, and experimental MaxMean each used 10,000 random gene sets per distinct pathway size, for short and middle phases, with score columns 5 and 6 analyzed separately in both languages. The earlier runs used 199 draws; these are new, full 10,000-draw runs.

## What ran

| Implementation | Score column | Configuration | Output directory |
|---|---:|---|---|
| R | 5 | `config.score5-full.R` | `results/r_score5_full` |
| R | 6 | `config.score6-full.R` | `results/r_score6_full` |
| Python | 5 | `config.score5-full.json` | `results/python_score5_full` |
| Python | 6 | `config.score6-full.json` | `results/python_score6_full` |

Each phase tests 2,638 short-term or 2,637 mid-term pathways. There are 372 distinct measured pathway sizes in short and 383 in middle. Sharing a null distribution among equal-size pathways gives **90,600,000 null-set evaluations**, **24 simulation result tables**, and **63,300 pathway-result rows** across the four configurations. Exact hypergeometric results are also saved for both phases in each configuration; this method does not need simulation. All four commands exited with status 0.

The runs use the staged primary CTgene lists, WITHOUTINF external data, and GO annotation. Membership counts are identical in Python and R: 39,958/40,235 measured genes and 29/138 measured CTgenes for short/middle. Score columns remain separate, larger supplied scores are assumed to indicate stronger evidence, and duplicate gene scores are explicitly averaged. These assumptions are recorded in each configuration.

## Verification

- All 17 staged/archive hashes match the input manifest.
- Every result table has the expected unique GO IDs, finite statistics, valid probabilities, and independently checked BH adjustment.
- Every simulated row records 10,000 draws and its seed. Its p-value matches `(null_exceedances + 1)/10001`.
- Python and R agree on scored genes, input summaries, pathway membership, counts, and deterministic statistics. The largest absolute observed-statistic difference was approximately `3.64e-12`. Exact hypergeometric probabilities agree within relative tolerance `1e-10`, including very small upper tails.
- Native Monte Carlo p-values differ, as expected from different random generators. The largest absolute difference was approximately 0.02770. No comparison fell outside the six-standard-error diagnostic screen; this is a broad diagnostic, not a formal equivalence test.
- An independent cross-language fixture enumerates all 120 three-gene subsets of ten genes and verifies all three statistics and the same finite-null p-value formula. Fixed and fitted FOSCO calculations agree too.
- All **12 Python tests** passed, including the complete synthetic analysis/correction/plotting commands. R numerical and synthetic integration checks passed. Both teaching examples ran successfully.
- Separate correction commands agree on the shared 20-gene synthetic input, both fitted lambdas, all 200 lambda-grid rows, and corrected scores. The largest table difference was approximately `1.34e-15`.
- Both plotting implementations produced eight figures per full configuration: **32 full-run figures**. Representative Python SUMSTAT/GSEA plots, the R MaxMean plot, and Python correction diagnostics were visually inspected.

Detailed numerical checks, per-method summaries, and SHA-256 fingerprints are in [full_simulation_verification.json](full_simulation_verification.json). Formula comparisons are in [cross_language_reference.json](cross_language_reference.json), and the correction command comparison is in [correction_language_comparison.json](correction_language_comparison.json). Four `*_full.log` files preserve the execution progress. Each result directory records its configuration and runtime environment.

## Fresh ZIP execution check

The updated ZIP was extracted to a new directory and **all 11 execution checks passed**: R numerical tests, R integration tests, the 12 Python tests, shared-null comparison, both teaching examples, both membership workflows, both plotting commands, and verification of every archived full-run table. The synthetic integration checks create new inputs and outputs. Full simulations were completed before packaging; the ZIP check verifies their saved results rather than repeating 90.6 million evaluations.

The checked code/configuration fingerprints, archive checkpoint hash, commands, exit codes, timings, and retained log location are in [latest_dual_language_package_check.json](latest_dual_language_package_check.json). The final ZIP adds this validation record and its report to the tested checkpoint; the tested code and configurations are unchanged.

## Numerical outcome snapshot

Counts below refer to **all eligible pathways with BH-adjusted p < 0.05**, before the plot's CTgene-overlap filter. They describe these runs and assumptions; they are not claims of replicated paper findings. Native Monte Carlo variability can change threshold membership.

| Score column | Phase | Method | R count | Python count |
|---:|---|---|---:|---:|
| 5 | Short | SUMSTAT | 0 | 0 |
| 5 | Short | GSEA | 0 | 0 |
| 5 | Short | MaxMean, experimental | 0 | 0 |
| 5 | Middle | SUMSTAT | 44 | 43 |
| 5 | Middle | GSEA | 39 | 39 |
| 5 | Middle | MaxMean, experimental | 40 | 41 |
| 6 | Short | SUMSTAT | 136 | 136 |
| 6 | Short | GSEA | 178 | 178 |
| 6 | Short | MaxMean, experimental | 157 | 157 |
| 6 | Middle | SUMSTAT | 208 | 217 |
| 6 | Middle | GSEA | 207 | 206 |
| 6 | Middle | MaxMean, experimental | 286 | 292 |

The exact hypergeometric counts are 11 short and 140 middle in both implementations and both score configurations. Changing a numerical score column does not change this membership-only test.

## Scope and environment

Tested locally on Windows with R 4.5.1, Python 3.11, NumPy 2.3.5, SciPy 1.17.1, and Matplotlib 3.10.7. R printed nonfatal startup locale warnings and the intended experimental-MaxMean notices; all commands completed. The saved R session uses locale `C`. GitHub Actions is configured but has not been run remotely.

The code now runs the requested short/middle statistical workflow in both languages. **Real-data FOSCO correction is still unconfirmed:** its raw score definition must be established before using the collaborator's values as correction inputs. Correction execution was tested on synthetic data. Score direction/construction and the scientific duplicate policy also need confirmation before interpreting the full simulations as final research results. This work does not reproduce every historical network, laboratory, supplementary-table, or exploratory analysis, and does not claim exact paper reproduction.

For readable explanations and runnable examples, start with [the walkthrough](../docs/WALKTHROUGH.md).
