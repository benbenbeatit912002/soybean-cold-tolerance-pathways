# Validation record — 15 October 2026

This is the **earlier R/199-draw validation record**. See [FULL_SIMULATION_VALIDATION.md](FULL_SIMULATION_VALIDATION.md) for the subsequent Python and R implementations, 10,000-draw runs, and cross-language checks. Statements below describe that earlier snapshot.

Environment: Windows, R 4.5.1 (ucrt); base R statistical functions only. Python audits use the standard library. GitHub Actions configuration is supplied but has not run remotely.

## Completed checks

| Check | Outcome |
|---|---|
| Parse every script in the original code archive | 9 parse; 2 fail. Exact failures saved in `legacy_parse.tsv`; parse success does not establish runnability |
| Parse all 10 new R/config/test files | Passed |
| Primary CTgene lists and GO annotation | 44 short / 143 middle / 170 union; 4,896 GO terms; all original declared term sizes match membership counts |
| External ID mapping | 1,388 older IDs resolved; 8 ambiguous and 7 blank-ID rows excluded with row-level records |
| Input provenance | 17 staged/archive hashes verified; originals also compared where present |
| Python import tests | Four tests passed: unique mapping, ambiguity, blank/unknown/unplaced IDs, and version vs capitalization handling |
| Numerical R checks | Passed: independent combinatorial hypergeometric tail, stable small-probability FOSCO, probability-domain guards, GSEA vs a full running sum, ties, duplicates, MaxMean, and seeded reproducibility |
| R integration checks | Passed: SUMSTAT, GSEA, and experimental MaxMean for both synthetic phases; separate both-phase correction entry point and diagnostic plot generation |
| Real-data hypergeometric analysis | Completed for both phases after audited mapping and deduplication of membership |
| Real-data score-column diagnostics | Completed both score columns × both phases × SUMSTAT/GSEA, at 199 null draws each, with explicit mean aggregation of duplicate scores |
| Output constraints | Unique GO rows, finite statistics, valid p/BH ranges, and Monte Carlo floor verified; see `output_verification.json` |
| Visual QA | Both hypergeometric plots and the synthetic middle correction diagnostic inspected; long pathway descriptions use explicit ellipses |

Resolved input and testing counts:

| Phase | Unique external genes | CTgenes measured / primary | Eligible pathways | Eligible pathways with CTgenes |
|---|---:|---:|---:|---:|
| Short | 39,958 | 29 / 44 | 2,638 | 96 |
| Middle | 40,235 | 138 / 143 | 2,637 | 335 |

These counts use the new 6–1,500 measured-gene size rule and full BH family. They are not the historical 100/342 candidate counts. Every missing CTgene is listed in the corresponding results TSV.

## Limits

The 199-draw score analyses are diagnostic only: the minimum nominal p-value is 0.005, and this resolution cannot support the paper's stringent adjusted-p selection. Do not interpret the absence of BH values below 0.05 in these runs as evidence that the methods or pathways are ineffective. A confirmed high-resolution rerun must follow decisions about score meaning, duplicate handling, and the scientific testing family.

The real-data FOSCO chain has been audited but not reconstructed from a confirmed raw-score definition. Its new entry point was validated on synthetic inputs. The original network, bar-chart group comparisons, wet-lab validation, and final supplementary result tables were not reproduced end to end. All 11 original code-folder files remain archived, including their known failures.

R initially inherited unsupported `C.UTF-8` locale settings and emitted startup warnings. The legacy parser sets a Windows UTF-8 character locale to enumerate Chinese filenames, and stops if no files are found. Clean analysis paths are ASCII. These warnings did not prevent the successful recorded runs.

Full run settings and `sessionInfo()` are saved beside each result. Generated outputs remain local in `results/`, which is ignored by Git. The local ZIP includes them for review, as well as this validation record. No upload or public release was performed.
