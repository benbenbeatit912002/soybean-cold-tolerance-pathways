# Short-score recalculation from upstream inputs — 21 October 2026

The Python and R implementations now **calculate scores from expression data,
gene correspondence, gene lengths, and SNP marker blocks**. Neither calculator
reads `gene_p_len.txt`, saved corrected probabilities, or saved short scores.
The old `--short-reference` score-copying option has been removed from both
correction commands. Its earlier preservation report is an archived record,
not evidence of recalculation.

**Numerical result:** the archived calculation plus reconstructed export rules
matches 42,574 of 42,575 saved short scores exactly at their stored precision.
An explicitly optional repair hypothesis also calculates the remaining score,
giving 42,575 exact decimal matches in both languages. **The historical origin
of that last repair is still unconfirmed.** Equal numbers do not establish that
the inferred repair was the procedure originally used.

## What was reconstructed

1. Read columns 5 and 8 of `耐寒資料庫0 1 24hr(轉換後)培修版.txt` as the short and
   middle input values. Preserve the source decimals. These are the supplied
   p-value fields; this branch does not calculate them from replicate counts.
2. Match either of the first two columns of `Glyma_genes_correspondence_final.txt`.
   Keep a record only when exactly one correspondence **row** matches and its
   target is nonempty. Two matching aliases in the same row are one match.
3. Join lengths from `Glyma2.0.txt`. Recount every `AX-...` marker between each
   gene header and `END` in `snps_allgenes_20kb.set.txt`.
4. Reproduce the archived 2,000-bp bins and their means **including zeros**.
   The archived breakpoints end at 98,001 and have no global fallback.
5. Search all 100 lambda values from 0.01 through 1.00, minimizing absolute
   Spearman correlation with SNP count. Keep duplicate records in their source
   order. The short-term optimum is 0.01.
6. Evaluate `1 - (1 - p)^(SNP^lambda)` using ordinary subtraction, then reproduce
   the saved probability precision and logged-score export described below.

The mapping retains 42,575 records representing 42,567 unique genes. It excludes
6,018 unmatched, 1,183 ambiguous, and one blank-target record from the 49,777
expression records. All reconstructed input fields and historical imputed SNP
counts agree with `gene_p_len.txt`. There are 2,569 zero counts before the
historical imputation and one afterward.

## Why the earlier recalculation differed

| Difference | Evidence and effect |
|---|---|
| Numerically stable FOSCO versus the archived expression | `log1p`/`expm1` avoids cancellation but changes the old tiny probabilities. It causes 55 additional short-score differences above `5e-8`. Ordinary subtraction resolves all 55. The largest is about 0.003742 for `Glyma.14G041700`. |
| Intermediate probability precision | The saved R table uses 15 significant digits. Re-importing that precision before the logarithm resolves a final rounding boundary for `Glyma.07G251000`. |
| Zero-score replacement | All 1,499 probabilities equal to one have saved logged scores of `1e-16`, rather than zero. This rule is observed in the saved artifacts; its original editing command was not found. |
| Final decimal export | The saved short scores use two decimals in scientific notation below `1e-4`, nine decimal places below 10 otherwise, and eight decimal places at or above 10, with trailing zeros omitted. This inferred rule matches every ordinary short score. |
| Different zero-SNP policy | The alternative nonzero group/global mean produces SNP 11 from the already-imputed mapped file. It does not reproduce the exceptional saved SNP 12 or its score. |

The inferred export rules are applied uniformly, without consulting any
reference score or its per-row number of digits. Both the unrounded logarithm
and the exported score are retained. The comparison uses **exact decimal
equality**, not a permissive error tolerance.

## The remaining gene and the average-value lead

`Glyma.14G129500` has length 52,247, raw short input `0.804571275`, and an empty
SNP block. It is the only mapped record in `[52001,54001)`, so the original
mapped-record imputation leaves it at zero. The older `1209` corrected file
also contains zero; a later corrected file instead records SNP 12 and short
probability `0.9947477`. The associated logged score is `0.002287056`.

The author recalls using an average, but did not confirm a population or the
lambda used for this row. Whole-genome length and SNP files contain two positive
donors in the same bin:

| Gene | Length (bp) | SNP markers |
|---|---:|---:|
| Glyma.06G180700 | 52,210 | 16 |
| Glyma.10G123400 | 52,443 | 7 |
| Glyma.14G129500 | 52,247 | 0 |

The positive-donor mean is `(16 + 7) / 2 = 11.5`, rounded to 12 by R and Python.
This is a plausible average-based explanation, **not a uniquely proven one**.
Other population/band choices can also yield 12; the mapped global positive
mean instead rounds to 11. Candidate means are recorded in the local audit
outputs.

SNP 12 alone does not resolve the score:

| Calculation | Short score before final decimal export |
|---|---:|
| SNP 12, fitted short lambda 0.01 | 0.09021082063909275 |
| SNP 12, fitted middle lambda 0.47 | 0.0022870490542424055 |
| Same middle-lambda calculation, probability rounded to 7 significant digits before the logarithm | 0.0022870563318260653 |

The last line exports to the original `0.002287056`. Its probability is computed
as `1-(1-0.804571275)^(12^0.47) = 0.9947477166692417`, then rounded to `0.9947477`.
Using the middle input in that same calculation also gives the saved middle
probability after 7-significant-digit rounding, providing a second numerical
check of the hypothesis.

The optional `--zero-repair-hypothesis` therefore calculates the whole-genome
positive bin mean for rows still at zero, fits the middle-phase lambda from the
data, and applies that lambda and 7-significant-digit rounding to those rows.
It does **not** contain a special gene ID, the constant SNP 12, a fixed lambda
0.47, or a copied probability/score in the calculation code. Ordinary short
records continue to use their fitted short lambda.

This is consistent with a manual repair after the middle-phase calculation
using the last lambda left in the R session and a printed console value. No
saved script or session proves that sequence. It is not a recommended new
statistical method and is not the default calculator behavior.

## Commands and verification

Run from the repository root. The four upstream research files must be present
under the specified original project directory.

```sh
python scripts/reproduce_historical_short.py "D:/Bioinformatics project" results/short_derivation_2026-10-21/python
Rscript --vanilla scripts/reproduce_historical_short.R "D:/Bioinformatics project" results/short_derivation_2026-10-21/r
```

Both write all outputs and return status 2 because the strict archived procedure
leaves one infinite score. PowerShell's command wrapper can report that as a
generic nonzero status. The incomplete row is listed in `unresolved_rows.tsv`.

To test the explicitly unconfirmed repair:

```sh
python scripts/reproduce_historical_short.py "D:/Bioinformatics project" results/short_derivation_2026-10-21/python_hypothesis --zero-repair-hypothesis
Rscript --vanilla scripts/reproduce_historical_short.R "D:/Bioinformatics project" results/short_derivation_2026-10-21/r_hypothesis --zero-repair-hypothesis
python scripts/verify_short_recalculation.py "D:/Bioinformatics project" results/short_derivation_2026-10-21
```

The verifier alone reads original expected results. It validates source hashes,
checks input reconstruction, aligns all duplicate records by occurrence, and
compares exported scores as decimals. See the
[machine-readable verification](short_recalculation_verification_2026-10-21.json).
The result report keeps `historical_zero_repair_provenance_confirmed: false`.

The new tests run both calculators in a synthetic directory containing only the
four upstream files. Adding deliberately invalid saved-score files leaves
results unchanged. Changing the exceptional raw short input changes its
calculated score. These tests distinguish calculation from reference copying.

Validation completed: all 23 Python tests passed, including integration tests
that execute both languages. The five historical-calculation tests were rerun
after the final input-serialization changes and passed. R core, integration,
reference-identity, and historical-calculation suites all passed. The final
independent real-data comparison checked four fresh runs (strict and hypothesis
in both languages). It found exact input reconstruction for all 42,575 rows and
exact decimal score equality for all rows in both hypothesis runs. Unrounded
corrected probabilities differ between languages by at most `1.11e-16`.

## Sources checked and remaining limits

The investigation checked the import/correction scripts, older saved result
versions, both saved R workspaces and available R history, original treatment
workbooks, and the draft corrected-score workbook. No repair command was found.
The draft corrected-score workbook contains cached values, with the exceptional
row still marked `#NUM!` and SNP zero. The workspaces contain older mapping or
later plotting objects rather than the required correction-session state.

The paper's [Methods](https://link.springer.com/article/10.1007/s00299-025-03643-2)
describe a fold-change-derived input, `10^(-abs(ln(FC)))`. The archived branch
reconstructed here reads supplied p-value columns. Numerical agreement with the
author-designated saved short scores does not reconcile that methodological
difference. This run does not regenerate the 10,000-draw pathway simulations or
establish reproduction of all published tables.

Original research files were not modified. The new files and reports are local;
this investigation did not publish or push a repository update.
