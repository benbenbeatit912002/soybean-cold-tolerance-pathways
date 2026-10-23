# Original short-score preservation — 21 October 2026

**Archived and superseded:** the current workflow recalculates scores from upstream inputs. The score-copying CLI option described below has been removed. This report records the earlier preservation check only; it is not the current workflow. See [actual recalculation and remaining provenance limits](SHORT_SCORE_RECALCULATION_2026-10-21.md).

**Result: PASS.** Both Python and R preserve all **42,575 original short score strings** and all **42,575 stored short adjusted-value strings** exactly. These archived short scores were used as the paper-result reference in this preservation check. The selected short fields now use that saved reference; the newly calculated values remain available in separate diagnostic fields.

This validates preservation of the author-designated original scores. It does not claim to regenerate every historical operation, reproduce all paper tables, or rerun the 10,000-draw pathway simulations. Machine-readable results, source hashes, and output hashes are in [paper_short_reference_verification_2026-10-21.json](paper_short_reference_verification_2026-10-21.json).

## Source and record identity

The reference contains the exact `p_1` score text from `Yamasaki_Go校正後1213_G_log10.txt`, together with the exact short adjusted value from `Yamasaki_Go校正後1213.txt`. The legacy SUMSTAT script reads the logged file at lines 3–4 and sums its second column at line 21. Hypergeometric, MaxMean, and one GSEA script also explicitly read this logged artifact.

The prepared input and reference are [`data/reference/historical_input.tsv`](../data/reference/historical_input.tsv) and [`data/reference/paper_short.tsv`](../data/reference/paper_short.tsv). Their companion `.provenance.json` files record source hashes and interpretation. Original research files were not changed.

All three historical files contain 42,575 records and 42,567 canonical gene IDs. Eight gene IDs occur twice, with conflicting scores but identical length and SNP count. A gene-only join, or even a gene/length/SNP join, would therefore be ambiguous. The raw and full-corrected files retain identical gene/length order; the logged file is reordered by decreasing length.

Reference staging joins by canonical gene ID, length, and occurrence order within each pair. Independent comparison with `-log10` of the stored corrected value validates every short match within `5e-8` (maximum `4.3841332514e-8`); the saved SNP counts also agree. This numerical tolerance validates correspondence, while the actual reference stores the original score strings without recalculation. All original records are retained one to one.

The correction commands then match reference rows by the unique physical `source_row` in `gene_p_len.txt`, requiring the same complete row set and checking gene ID, raw short value, length, and original SNP count. `reference_log_row` records the corresponding physical line in the logged source. Staging checks pinned source hashes and refuses to replace an existing different reference.

## Selected scores and new calculations

| Output field | Meaning with `--short-reference` |
|---|---|
| `score_short` | Original logged score, copied as exact stored text |
| `adjusted_short` | Original adjusted short value, copied as exact stored text |
| `recomputed_score_short`, `recomputed_adjusted_short` | Independent new calculation after the alternative SNP imputation |
| `snp_original`, `snp` | Saved input count and current imputed count |
| `reference_snp_short` | Count stored with the historical short reference |
| `short_score_source` | `historical_reference` |
| `score_middle`, `adjusted_middle` | New middle-phase calculation after imputation |

The fitted lambda and search grids describe the recomputation: **0.01 short and 0.47 middle**. The selected historical short values are not claimed to follow that new fit. In `correction_summary.tsv`, `adjusted_spearman` uses selected adjusted values against current imputed counts, while `recomputed_adjusted_spearman` uses the newly calculated values. Short diagnostic plots explicitly show the recomputation.

For **`Glyma.14G129500`**, source line 29,915 maps to logged line 11:

| Quantity | Value |
|---|---:|
| Saved input SNP | 0 |
| New imputed SNP | 11 |
| Historical reference SNP | 12 |
| Preserved adjusted short value | `0.9947477` |
| Preserved short score | **`0.002287056`** |
| Recomputed adjusted short value | `0.8121628169781555` |
| Recomputed short score | `0.09035689757297066` |

The new SNP count still uses the same-group positive mean, with the overall positive mean as fallback. This gene has no positive donor in its length group, so the overall mean `10.731667214731996` rounds to 11. That rule remains available for reanalysis; it does not overwrite the selected historical short score.

**Numerical clue, not recovered provenance:** applying lambda **0.47** to this gene's raw short value and SNP **12** gives `0.9947477166692417`, which rounds to the stored adjusted value `0.9947477`. This is consistent with the shape of the exceptional value, but no recovered command establishes that this operation generated it. Re-logging that full-precision model value would also differ from the saved logged string. The reference therefore preserves the original string rather than substituting an inferred formula.

## Verification results

| Check | Result |
|---|---:|
| Records retained in each language | 42,575 |
| Unique canonical genes | 42,567 |
| Exact original short score text matches per language | **42,575 / 42,575** |
| Exact original adjusted short text matches per language | **42,575 / 42,575** |
| Maximum Python/R selected short score difference | **0** |
| Maximum Python/R selected adjusted short difference | **0** |
| Maximum Python/R recomputed short score difference | `4.98e-14` |
| Maximum Python/R recomputed middle score difference | `4.98e-14` |
| Recomputed versus saved short score differences greater than `5e-8` | 56 |

The new calculation differs from saved short scores through stable numerical arithmetic as well as the exceptional gene's historical handling. Preserving only that one gene would not preserve the entire original score vector. Both implementations therefore preserve every reference score. `short_reference_comparison.tsv` exposes all recalculated differences; their maximum absolute magnitude is `0.08806984157297067`.

Both real-data commands completed, and the independent verifier compared outputs with the original sources. All **18 Python tests** passed, as did the R numerical, integration, and reference tests. Reference tests cover selection and record-identity checks. Both short diagnostic images were visually checked: their headers clearly identify the recomputation and retained original scores, with readable labels and no clipping. Exact equality of selected scores is expected because both languages preserve the same source data; it is not independent evidence that the historical derivation was correct.

## Reproduce the check

Run from the repository root using the included prepared files:

```sh
python scripts/correct_scores.py data/reference/historical_input.tsv results/paper_short_reference_2026-10-21/python --duplicate-policy preserve --short-reference data/reference/paper_short.tsv
Rscript --vanilla scripts/correct_scores.R data/reference/historical_input.tsv results/paper_short_reference_2026-10-21/r --duplicate-policy preserve --short-reference data/reference/paper_short.tsv
python scripts/verify_paper_short_reference.py "D:/Bioinformatics project" results/paper_short_reference_2026-10-21
```

The verifier requires the original source files under its first argument. To restage the included historical input and reference from those originals:

```sh
python scripts/prepare_historical_correction.py "D:/Bioinformatics project" data/reference/historical_input.tsv --short-reference data/reference/paper_short.tsv
```

On the author's Windows installation, invoke R as `& 'C:/Program Files/R/R-4.5.1/bin/Rscript.exe'` in PowerShell if `Rscript` is not on `PATH`.

The [20 October correction run](SNP_IMPUTATION_VALIDATION_2026-10-20.md) remains a documented reanalysis. Its recomputed short scores are not the selected original scores. The existing WITHOUTINF simulations use a different score branch and do not read the reference automatically.
