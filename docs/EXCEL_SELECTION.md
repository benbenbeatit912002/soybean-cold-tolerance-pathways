# Programmatic selection from the Excel candidate tables

Excel conditions can be implemented as reproducible code. The Python and base-R selectors read complete candidate tables, calculate each condition, and export every inclusion/exclusion decision. They never read the expected selected list during calculation.

The packaged candidate table has **1,326 method rows**: three methods for 100 short-term and 342 mid-term GO terms. Each normalized row was checked against its exact source range in the local June 2025 `TableS2_Pathway enrichment analysis.xlsx`; all fields agree. It is a historical candidate table, not freshly recalculated enrichment output.

## Run either language

```sh
python scripts/select_pathways.py data/reference/crosstalk/selection_candidates.tsv config.selection-candidate.json results/selection_python
```

```sh
Rscript --vanilla scripts/select_pathways.R data/reference/crosstalk/selection_candidates.tsv config.selection-candidate.R results/selection_r
```

Choose a new output directory for each run. Python uses only the standard library for this step; R uses base R. Excel and its original drive paths are not needed.

Both included configurations implement the following **unconfirmed candidate reconstruction**:

| Setting | Short | Middle |
|---|---|---|
| Count source | Hypergeometric candidate table | Hypergeometric candidate table |
| Biological condition | CTgene count > 5 **OR** CTgene ratio >= 0.05 | CTgene count >= 40 |
| Statistical condition | At least one main method has nominal p < 0.05 | At least one main method has nominal p < 0.05 |
| Combination | Biological AND statistical condition | Biological AND statistical condition |

Counts are integers, so `min_ct: 6` implements `count > 5`. The ratio denominator is the annotated pathway size. The three methods are hypergeometric, SUMSTAT and GSEA; MaxMean does not participate.

These rules reproduce the archived **13 short / 11 middle IDs exactly**, when applied to the complete historical candidate table. The author could not confirm that these were the actual Excel operations, so this remains a candidate reconstruction. Matching a selected list does not uniquely identify the original rules, and the nominal-p threshold differs from a universal BH threshold described in the paper. The primary reference crosstalk workflow continues to use the preserved historical Excel list.

Changing short `> 5` to `>= 5` selects three additional terms: GO:0002679, GO:0035556 and GO:0043565. Substituting BH-adjusted values also changes the selection. These are meaningful choices, not formatting differences. The historical method tables disagree on some CTgene counts; the explicit `count_method` avoids silently taking whichever row is encountered first.

## Read or change the rules

| Configuration field | Meaning |
|---|---|
| `count_method` | Which method table supplies CT count and description |
| `probability_field` | `p_value` for nominal values, or `p_bh` for BH-adjusted values |
| `alpha` | Strict upper probability threshold: p < alpha |
| `method_rule` | `any` or `all` of the three methods |
| Per-phase `min_ct` | Inclusive minimum integer count |
| Per-phase `min_ratio` | Inclusive minimum CT count / annotation size |
| Per-phase `biological_rule` | `count`, `ratio`, `and`, or `or` |

Outputs are `selection_audit.tsv` for all candidates and `selected_pathways.tsv` for selected rows. Each row records counts, ratio, the actual probabilities used, condition flags, agreement between method counts and the selection label. Saved settings and input checksums preserve the configuration used.

For other data, supply a UTF-8 tab-separated table with `phase`, `method`, `go_id`, `description`, `annotation_size`, `ct_count`, `p_value` and `p_bh`. There must be one row for each phase/GO/method combination and all three methods per GO term. Duplicate rows, missing methods, inconsistent sizes and invalid probabilities are rejected. Use a descriptive label for your rules; do not label a new selection as the historical procedure.

The resulting phase/GO table can be supplied to either standalone crosstalk command using `--selected-pathways` and an explicit `--selection-label`. Its provenance remains separate from the freshly calculated enrichment probabilities displayed in crosstalk audits.

## Count lineage after selection

The selected IDs match the final original Excel, but two short-term count fields differ: GO:0045893 has 15 versus 14 CTgenes, and GO:0050832 has 12 versus 11. The candidate table counts the complete primary CT list; the final Excel counts the measured subset. Glyma.01G216200 is annotated to both pathways but is absent from the mapped expression input. Both main analysis implementations already exclude it from the measured CT foreground.

The [membership-lineage audit](../reports/CT_MEMBERSHIP_LINEAGE_2026-10-30.md) reproduces all 100 archived January short hypergeometric rows in Python and R and documents the remaining mixed manuscript versions. It also records 17 middle pathways with differing SUMSTAT/GSEA candidate counts. These findings do not confirm the secondary-selection rule and do not authorize replacing the archived candidate table with recalculated statistics.
