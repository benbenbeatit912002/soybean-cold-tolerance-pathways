# Fresh ZIP execution test — 15 October 2026

This record covers the earlier R/199-draw archive. The later Python/R full-run package is documented in [FULL_SIMULATION_VALIDATION.md](FULL_SIMULATION_VALIDATION.md).

**Result: all nine execution checks passed. No code fixes were needed.**

The delivered ZIP was extracted to a new folder. Existing results and the old output-verification report were excluded before running, so success required fresh output generation. The analysis used only the extracted package inputs.

| Check | Exit code | Seconds |
|---|---:|---:|
| core_checks | 0 | 0.29 |
| integration_checks | 0 | 0.73 |
| import_checks | 0 | 0.18 |
| legacy_parse_audit | 0 | 0.26 |
| hypergeometric_both_phases | 0 | 14.35 |
| score5_both_phases | 0 | 46.86 |
| score6_both_phases | 0 | 46.1 |
| plot_generation | 0 | 0.3 |
| output_validation | 0 | 0.16 |

All 33 regenerated result TSV files exactly match their archived counterparts by SHA-256. All 4 generated PNGs have valid signatures and nonzero dimensions.

Real-data coverage: short and middle hypergeometric analysis; supplied score columns 5 and 6 separately, each with SUMSTAT and GSEA at 199 random gene sets. Correction and experimental MaxMean were exercised with synthetic inputs by the integration test.

The legacy parse audit completed successfully as an audit: it still reports the two already documented original scripts with syntax errors. The runnable workflow does not execute those archived scripts. R emitted inherited C.UTF-8 startup-locale warnings; they did not prevent execution. The experimental MaxMean warning is intentional.

This verifies successful execution and repeatability on Windows with R 4.5.1. It does not establish exact scientific reproduction of the paper or validate the unresolved raw-score definitions.

Archive tested (SHA-256): `aad818a9c5ad238a535ca8030aa03eb78a3668300ee03560a522e701b6fd00b2`.

Raw logs and fresh outputs: `D:\Bioinformatics project\package-run-check-pbfuj0w5`.

Machine-readable record: [latest_execution_check.json](latest_execution_check.json).
