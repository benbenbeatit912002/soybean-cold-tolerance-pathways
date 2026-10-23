# Inputs and provenance

The current portable workflow is described in [data/README.md](../data/README.md). All required inputs are included in `data/upstream/` and `data/input/`; users do not need the original project directory. The source names below document the author's archive. Files under `data/input/` are prepared copies; the originals were not edited.

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

`external_*.tsv`: `gene_id`, `source_col_2`, ..., `source_col_6`. The original WITHOUTINF files have no headers. Columns 5 and 6 were supplied by Edward Kao. The 18 October 2026 provenance audit recovered historical headers and matched every raw WITHOUTINF numeric tuple to the original treatment workbook: columns 2–4 are treatment replicates and column 5 is fold change (treatment-replicate total / control-replicate total). Column 6 is historically labelled `Pvalue`, but exceeds one and its construction remains unresolved; it is not a raw probability. The historical score used by each method still needs to be established. Numeric columns are retained unchanged. WITHOUTINF denotes exclusion of records with no information; the audit independently verifies finite numeric values rather than inferring behavior from the folder name.

Header sources are `研究法/大豆基因/數據/外部資料匯入/舊版/耐寒資料庫 short _pvalue.txt` and `耐寒資料庫 middle _pvalue.txt`, line 1. The original workbook is `研究法/大豆基因/數據/外部資料匯入/Pvalue轉換_Yamasaki (2016) 耐寒資料庫0 1 24hr.xlsx`, sheets `1hr cold treatment` and `24hr cold treatment`, headers `A1:F1`. See the [paper concordance audit](../reports/PAPER_CONCORDANCE_2026-10-18.md) for matching tolerances and remaining differences.

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

The follow-up [data-flow audit](../reports/HISTORICAL_DATA_FLOW_2026-10-18.md) exactly replays the external mapping: 49,777 source rows become 42,576 after unmatched/ambiguous-row exclusions, then 42,575 after excluding one empty target. The saved full file's actual GO intersection is 37,153 unique genes. Both saved FOSCO phases are numerically reproduced for the 42,574 positive-SNP rows using lambda 0.01/0.47 and the source p-value columns. The Peixiu file also contains separate cached FC columns consistent with `10^logFC`; these are not the columns read by the original correction block.

The legacy `p_1` and `p_24` headers are reused for raw probabilities and transformed scores in different files. The logged file contains values greater than one; its columns must not be fed back through `-log10` as if they were probabilities. The report compares saved files without treating a matching filename as proof of matching semantics.

Hashes and transformations: `reports/source_manifest.json`. Raw input counts: `reports/input_audit.json`. Resolved import counts: `reports/import_summary.json`. Differences among saved correction artifacts: `reports/lineage_audit.json`.

## Correction input and imputation audit

Both correction entry points read the following tab-separated columns:

| Column | Meaning |
|---|---|
| `gene_id` | Gene identifier; normalized to the package's canonical capitalization |
| `raw_short`, `raw_middle` | Values in `(0,1]` before FOSCO and before `-log10` |
| `snp` | Finite nonnegative integer count; zero invokes the documented imputation rule |
| `gene_length` | Positive integer length in bp, required when any SNP count is zero |
| `source_row` | Original physical line number, preserved in outputs; used when independently comparing historical records |

`scripts/prepare_historical_correction.py` stages the saved `研究法/大豆基因/數據/外部基因校正/gene_p_len.txt`, mapping `p_1`, `p_24`, and `len` to the named columns above. It preserves all 42,575 records and writes a source hash. **This saved file already contains earlier SNP imputations**; staging it does not recover the original marker counts or recompute those previous replacements.

Correction outputs retain `snp_original`, put the final count in `snp`, and label each row with `snp_imputation_source`. `snp_imputation.tsv` reports data-row number (excluding the header), gene ID, bin bounds, original/replacement counts, positive donor count, unrounded mean, and whether the group or global mean was used. The 20 October rerun preserves duplicate records explicitly for comparison. See the [validation report](../reports/SNP_IMPUTATION_VALIDATION_2026-10-20.md) for results and limits.

## Original short-score reference

The original short scores serve as the archived comparison reference for the paper's short-term analysis. The prepared files `data/reference/historical_input.tsv` and `data/reference/paper_short.tsv` retain all 42,575 historical records, including the eight extra duplicate rows. Each has a companion `.provenance.json` file recording its source and scope. These files belong to the historical correction branch, separate from `data/input/external_*.tsv` and the WITHOUTINF simulations.

`paper_short.tsv` contains these columns:

| Column | Meaning |
|---|---|
| `source_row` | Unique physical line in `gene_p_len.txt`; binds the reference to one original record |
| `gene_id`, `raw_short`, `gene_length`, `snp` | Original identity and input values, checked before the reference is used; `snp` is the count before the current imputation |
| `reference_adjusted_short` | Exact saved `p_1` text from `Yamasaki_Go校正後1213.txt` |
| `reference_score_short` | Exact saved `p_1` score text from `Yamasaki_Go校正後1213_G_log10.txt` |
| `reference_snp_short` | SNP count stored with the logged reference, which can differ from the original and currently imputed counts |
| `reference_log_row` | Physical line in the reordered `_G_log10` file |

Staging uses canonical gene ID, length, and occurrence order within each pair to align the logged file, then independently validates the score and SNP correspondence. A gene-only join, or even a gene/length/SNP join, cannot distinguish the eight duplicate pairs. The source hashes prevent silently rebuilding this reference from changed historical files.

The references are expected results for independent verification only. The score-copying option has been removed. `reproduce_historical_short.py` and `.R` do not read either prepared reference file; they reconstruct the raw inputs and calculate from the four original upstream files.

Historical recalculation output fields include `source_row`, `expression_row`, `original_gene_id`, `gene_id`, `raw_short`, `raw_middle`, `gene_length`, `snp_markers`, `historical_snp`, `snp`, `adjusted_short`, `adjusted_short_saved_precision`, `score_short_unrounded`, `score_short`, and `calculation_status`. Raw input decimals are preserved as supplied. `snp_markers` is recounted from raw blocks. `historical_snp` is the old mean-including-zeros result. `snp` differs only if the optional repair hypothesis is enabled. The unrounded score is the logarithm after probability serialization and the observed zero-score replacement, before final decimal formatting.

The strict procedure leaves `Glyma.14G129500` at zero SNP. The optional hypothesis derives SNP 12 from whole-genome donors 16 and 7 and calculates its saved short score using the fitted middle lambda and seven-significant-digit rounding. Its provenance remains unconfirmed. See the [recalculation report](../reports/SHORT_SCORE_RECALCULATION_2026-10-21.md).

## Fresh two-phase analysis outputs

`recalculated_scores.tsv` preserves all 42,575 records. `score_short` and
`score_middle` are the finite, exported statistics consumed directly by
`config.recalculated-full.json` / `.R`. They are already transformed scores;
do not take another logarithm. The corresponding `adjusted_*` fields preserve
the calculated probabilities, including a cancellation-induced zero whose
logged score requires an explicitly labeled repair. The `score_*_unrounded`
field contains the logarithm before normal export rounding; for the exceptional
infinite-score replacement it contains the repaired finite score.

`snp_markers`, `historical_snp`, and `snp` distinguish original marker counts, the
old bin-mean imputation, and the optional exceptional repair. `status_short`
and `status_middle` identify calculated, hypothesized, and unresolved records.
`exception_audit.tsv` records the derived maximum finite score and every
exception. These fields do not confirm that the hypothesized manual repairs
were performed historically.

The new analysis configuration uses `score_input` and `score_field`, with
`{phase}` replaced by `short` or `middle`. It averages duplicate records only
when constructing the unique-gene analysis universe. The additional historical
SUMSTAT replay keeps those records separate. Full settings and comparison
limits are in [the rerun report](../reports/RECALCULATED_ANALYSIS_2026-10-21.md).
