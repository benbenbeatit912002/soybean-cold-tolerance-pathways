# Historical Excel pathway selection

These inputs preserve the selected pathway identities and provide independent checks. They do not supply gene scores or replace current enrichment results.

- `selected_pathways.tsv`: 13 short and 11 middle selections, original counts/BH values, sheet and row.
- `selection_criteria_audit.tsv`: separate count, ratio and BH checks on historical values; no inferred filter is applied.
- `saved_jaccard_pairs.tsv`: 78 short and 55 middle unordered pairs from archived February 21 matrices, used only for validation.
- `provenance.json`: workbook ranges, original-file SHA-256 hashes and checked manuscript/list sources.

Original workbooks are not runtime dependencies. See [crosstalk documentation](../../../docs/CROSSTALK.md) for interpretation and commands.
