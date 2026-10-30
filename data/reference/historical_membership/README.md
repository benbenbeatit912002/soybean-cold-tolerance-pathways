# Historical membership validation references

These are validation-only snapshots of the archived January 09 hypergeometric tables and the short-term scored CT list. The Python and R membership audit commands do not read these files during calculation. They derive membership from the primary CT lists, the GO annotation and the supplied gene universe.

- `short_0109.tsv`: 100 original short-term rows, including the saved probability strings and their precision.
- `middle_0109.tsv`: 342 original mid-term rows.
- `short_scored_genes.tsv`: the 41 IDs in the archived scored CT list; numerical gene scores are deliberately omitted.
- `provenance.json`: original relative source paths, checksums and normalization details.

The original tables use full annotated pathway sizes, a measured gene universe, a measured CT foreground and `1 - phyper(...)`. This is a historical convention, not the main analysis's consistent measured-membership test. The stable upper tail is also exported by the audit. Near machine epsilon, subtraction results can differ between numerical libraries even when their stable tails agree.

See [the count-lineage report](../../../reports/CT_MEMBERSHIP_LINEAGE_2026-10-30.md) for results, unresolved table differences and commands.
