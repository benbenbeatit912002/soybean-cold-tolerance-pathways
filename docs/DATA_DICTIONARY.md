# Inputs and provenance

All source paths below are relative to the original `D:/Bioinformatics project` directory. Files under `data/input/` are prepared copies; the originals were not edited.

| Source | Staged file | Meaning |
|---|---|---|
| `研究法/大豆基因/CTgenes 含分數/SHORT.txt` | `ct_short.tsv` | Primary short-term CTgene membership and prior score |
| `研究法/大豆基因/CTgenes 含分數/MIDDLE.txt` | `ct_middle.tsv` | Primary mid-term CTgene membership and prior score |
| `研究法/大豆基因/大豆基因 外部/WITHOUTINF/outer_data_2.0_short_address.txt` | `external_short.tsv` | External short-term rows after audited ID resolution |
| Same directory, `outer_data_2.0_middle_address.txt` | `external_middle.tsv` | External mid-term rows after audited ID resolution |
| `研究法/大豆基因/數據/GO_terms_4899.txt` | `pathways.gmt` | GO IDs, descriptions, and gene membership |
| `研究法/大豆基因/大豆基因 外部/Glyma_genes_correspondence_final.txt` | `gene_correspondence.tsv` | Wm82.a1.v1 / a1.v1.1 to a2.v1 mapping |

## Column definitions

`ct_*.tsv`: `gene_id`, `prior_score`. The supplied primary lists have no duplicate IDs. The short-term score range is 8–12.2; mid-term is 10.03511753–17.80937504. These are not probability values. This package preserves them without inventing a transformation.

`external_*.tsv`: `gene_id`, `source_col_2`, ..., `source_col_6`. The original files have no headers. The author confirmed that columns 5 and 6 are scores supplied by Edward Kao. Units, score formulas, directionality, meanings of columns 2–4, and which score entered each historical method still need documentation. Numeric columns are retained unchanged. The author described WITHOUTINF as excluding data with no information; the audit independently verifies finite numeric values rather than inferring behavior from the folder name.

`pathways.gmt`: one tab-separated record per GO term: GO ID, description, then all member gene IDs. The original `N.genes` values agree with the count of unique member IDs in every row. The table contains 4,896 terms and 44,964 unique annotated genes; it does not itself contain per-gene numerical scores.

## Import checks

| Check | Short | Middle |
|---|---:|---:|
| Raw external rows | 40,061 | 40,332 |
| Old-format IDs uniquely remapped | 1,388 | 0 |
| Ambiguous mappings excluded | 8 | 0 |
| Blank-ID rows excluded | 5 | 2 |
| Retained external rows | 40,048 | 40,330 |
| Unique resolved gene IDs | 39,958 | 40,235 |
| Extra duplicate rows after resolution | 90 | 95 |

Unplaced `Glyma.u......` IDs are retained when already in the supplied v2 format. Capitalization of the chromosome-associated `g/G` is normalized. This normalization alone is not a genome-version conversion; older identifiers require the correspondence table. No ambiguous target is selected arbitrarily.

`reports/*_import_decisions.tsv` gives original row numbers, original IDs, resolved IDs, and outcomes. Duplicate score rows are retained in staged inputs so that a later aggregation policy remains visible. Membership analysis uses each gene once. Diagnostic score runs explicitly take the arithmetic mean of the supplied score for each resolved gene; the main configuration rejects duplicates for score analyses until a policy is chosen.

## Historical files are a different data branch

The legacy statistics scripts load `Yamasaki_Go校正後1213_G_log10.txt`, which has 42,575 rows but only 42,567 unique IDs, including eight conflicting duplicated genes. This is not either named WITHOUTINF input. `short_CT_G_score.txt` contains 41 genes, whereas the primary short list contains 44. The three missing IDs are listed in `reports/input_audit.json`.

The legacy `p_1` and `p_24` headers are reused for raw probabilities and transformed scores in different files. The logged file contains values greater than one; its columns must not be fed back through `-log10` as if they were probabilities. The report compares saved files without treating a matching filename as proof of matching semantics.

Hashes and transformations: `reports/source_manifest.json`. Raw input counts: `reports/input_audit.json`. Resolved import counts: `reports/import_summary.json`. Differences among saved correction artifacts: `reports/lineage_audit.json`.
