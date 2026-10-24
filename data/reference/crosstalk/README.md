# Historical Excel pathway selection

These inputs preserve the selected pathway identities and provide independent checks. They do not supply gene scores or replace current enrichment results.

- `selected_pathways.tsv`: 13 short and 11 middle selections, original counts/BH values, sheet and row.
- `selection_criteria_audit.tsv`: separate count, ratio and BH checks on historical values; no inferred filter is applied.
- `saved_jaccard_pairs.tsv`: 78 short and 55 middle unordered pairs from archived February 21 matrices, used only for validation.
- `provenance.json`: workbook ranges, original-file SHA-256 hashes and checked manuscript/list sources.
- `trace_provenance.json`: dated original network script, cutoffs, precision and gene-group source hashes.
- `historical_middle_14.tsv`: GO sets identified by exact gene membership in the earlier middle export.
- `saved_middle_122_ctgenes.tsv`: normalized archived middle CT list, used only for checking recalculated membership.
- `selection_candidates.tsv`: all three main methods for all 100 short and 342 middle candidates, with source sheet/range per row.

Original workbooks are not runtime dependencies. See [crosstalk documentation](../../../docs/CROSSTALK.md) for interpretation and commands.
