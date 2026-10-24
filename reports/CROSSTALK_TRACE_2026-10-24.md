# Resolving edge counts, middle CT membership and Excel selection

This follow-up supersedes the unresolved edge-count explanation in the earlier [crosstalk report](CROSSTALK_2026-10-24.md). It recovers an additional dated original script and traces gene membership through the saved pathway groups. Original research files were read without modification.

## Why 67 edges became 61

The earlier reanalysis interpreted the paper's fractions as retaining the largest 80% of positive short pairs and 50% of positive middle pairs, rounding upward and keeping ties. That operational definition retains 39 and 28 edges, or 67 after merging.

The original `crosstalk 20240916.R` instead reads the final February 21 matrices and sets values below fixed cutoffs to zero: **0.004456328** for short (line 31), **0.14160401** for middle (line 73). The older undated network script reads an earlier 14-pathway middle matrix with a different cutoff; it is not the correct final script.

The matrices contain Jaccard values printed to nine decimal places. Applying the dated script's cutoffs at that precision retains **36 short + 25 middle = 61 merged edges**. Both current implementations now recalculate intersections/unions from annotation, round only the value used for filtering, and apply those source cutoffs. No desired edge count is an algorithm input.

Rounding affects a boundary edge: GO:0009631–GO:0045893 has J = 5/1122 = 0.004456327985739751. Comparing full precision with 0.004456328 would exclude it; comparing the saved nine-decimal value includes it. The raw coefficient and `filter_jaccard` are both exported.

The two cutoffs equal the nine-decimal medians of their complete square matrices, including diagonal ones, and the dated script comments indicate 0.5. This explains a compatible origin of the constants. A complete saved quantile-calculation expression was not found; the verified replay uses the explicit constants rather than claiming an unobserved formula.

The six edges retained by the fraction interpretation but excluded by the dated script are:

| Phase | Source | Target | Jaccard |
|---|---|---|---:|
| Short | GO:0009737 | GO:0035264 | 0.004003203 |
| Short | GO:0010117 | GO:0010200 | 0.004393673 |
| Short | GO:0010200 | GO:0035264 | 0.004424779 |
| Middle | GO:0009414 | GO:0009867 | 0.130551817 |
| Middle | GO:0009620 | GO:0009867 | 0.129591837 |
| Middle | GO:0009723 | GO:0010363 | 0.127839066 |

## The additional historical middle CTgene

The difference is **Glyma.09G149200**. It is present in the primary middle CT list, so it was not lost during CT import or score reconstruction.

The saved `cluster2.txt` contains **14 pathway gene sets**. Each row matches a complete packaged GO set exactly. The final middle selection contains only 11 of these. The removed pathways are GO:0009793 (embryo development ending in seed dormancy), GO:0009826 (unidimensional cell growth), and GO:0048046 (apoplast).

Glyma.09G149200 occurs in the first two removed pathways and in none of the final 11 middle pathways. Recomputing the union of the 14 saved pathway sets and intersecting with the same primary middle CT list yields exactly the **122 IDs** in `ct_cluster2_in_list.txt`. Repeating the operation with the final 11 sets gives 121, with Glyma.09G149200 as the only difference.

This establishes a pathway-version mismatch between the saved middle gene export and final pathway list. The final 11-pathway output remains 121; the gene is not inserted artificially. It remains in the combined export because a selected short pathway, GO:0009737, contains it and combined CT membership uses both phase CT lists. Thus the combined total can remain 145.

## Excel rules

Both language implementations can now select from all 100/342 candidate pathways, with an audit of every condition. A candidate rule reproduces 13/11 exactly: short `(count > 5 OR ratio >= 0.05)` and middle `count >= 40`, each combined with any main method's nominal p < 0.05. Counts come explicitly from the hypergeometric candidate sheet. The author requested that this remain unconfirmed. It is not silently promoted to the original method or used to replace the reference list.

All 1,326 normalized method rows were checked against the source workbook cells with zero differences. Python and R produce the same selection audits. See [selection commands and definitions](../docs/EXCEL_SELECTION.md).

## Repeat the checks

```sh
python scripts/trace_crosstalk.py results/crosstalk_lineage
```

This independently recomputes both edge procedures, exports their six-edge difference, and checks all historical middle CT identities against the archived list. The evidence is packaged under `data/reference/crosstalk/`, with hashes in `trace_provenance.json`.

Both standalone network runs completed on the existing full enrichment results. Nine deterministic tables and eight gene lists agreed across languages. All 133 recomputed Jaccard pairs match archived precision; the retained sets also exactly match application of the dated script to the saved matrices. Unit checks exercise the rounding boundary, selection inequalities and input validation. No enrichment simulations or external SoyNet queries were required for these downstream calculations.
