# Original script data dependencies — 18 October 2026

This inventory covers every file found recursively in `D:/Bioinformatics project/研究法/大豆基因/程式`: 11 R scripts stored as `.txt` or `.R`. All 11 are byte-for-byte identical to their copies in the GitHub package’s `legacy/` directory; no original script in this folder was omitted. It does not claim that code elsewhere in the project has been inventoried.

**Status: the declared input/output calls are mapped; complete upstream data lineage is not yet resolved.** These are static code observations, not execution results. The original scripts were not executed and source/data files were not changed. R `parse()` was used only to check syntax.

The scan found 39 active `read.table()` calls, 15 active `write.table()` calls, one commented-out write, no `source()` calls, and no explicit plot-file export calls. Counts describe call sites, not distinct files, successful reads, or successful outputs.

`L4` means source line 4. The table shows basenames for readability. Exact original paths, line contents, inferred current paths, and alternative basename matches are recorded in:

- [Full machine-readable inventory](script_io_inventory_2026-10-18.json)
- [Every input/output call](script_io_inventory_2026-10-18.tsv)
- [Per-script summary](script_summary_2026-10-18.tsv)
- [Parse-only results](original_script_parse_only.tsv)

Original paths mostly point to `C:/Users/chris/Desktop/程式語言/研究法/大豆基因`. For existence checks, this prefix was replaced with the current project folder. This is an audit inference; it has not repaired the original code. A file existing at the inferred path does not prove it is the version used for publication.

## All 11 scripts

| Script | Declared file inputs | Declared outputs / plots | Session dependencies and gaps |
|---|---|---|---|
| `#外部資料匯入與校正繪圖.txt` | `耐寒資料庫0 1 24hr(轉換後)培修版.txt` (L4); `Glyma_genes_correspondence_final.txt` (L8); `Glyma2.0.txt` (L12); `snps_allgenes_20kb.set.txt` (L17); `GO_terms_4899.txt` (L81); `answer1.txt` (L97); `answer2.txt` (L99) | `gene_p_len.txt` (L75); `gene_p_len_go_yamasaki.txt` (L92); `Yamasaki_Go校正後1213.txt` (L174) Creates plot1 at 140, but contains no explicit file export. | No undeclared input object identified in this static pass. |
| `crosstalk 20240916.R` | `jc_point_short.txt` (L8); `jc_point_middle.txt` (L12); `both.txt` (L15); `only_short.txt` (L19); `only_middle.txt` (L21) | No table file written. Opens tkplot at 163; contains no explicit saved plot path. | No undeclared input object identified in this static pass. |
| `gsea_菲斯福氣.txt` | `short_CT_G.txt` (L27); `answer2.txt` (L31); `short_CT_G_score.txt` (L36) **path unresolved**; `middle_CT_G_score.txt` (L38) **path unresolved**; `Yamasaki_Go校正後1213_G_log10.txt` (L42) | `middle_esmax.txt` (L87) | `dat` (L3); `CTgenes_short` (L16): Used before it is loaded at 27.; `all_genes_short_outside_log` (L157); `number_sample` (L219): Subassignment before initialization in this file. |
| `gsea全部.txt` | `answer1.txt` (L3) **path unresolved**; `answer1.txt` (L49) **path unresolved** | `short_esmax.txt` (L45) | `all_genes_short_outside_log` (L5); `sample_intersect_2` (L88): Used before its later assignment at 97.; `和` (L91): Bare symbol, not a comment.; `o` (L118): Bare symbol, not a comment. |
| `maxmean0127版本.txt` | `Yamasaki_Go校正後1213_G_log10.txt` (L7); `answer1.txt` (L11); `answer2.txt` (L83) | `maxmean_short_p0127.txt` (L70); `maxmean_middle_p0127.txt` (L133) | `k1; intersection_fit_data` (L66); `k2; intersection_fit_data2` (L129) |
| `maxmean原版.txt` | `Yamasaki_Go校正後1213_G_log10.txt` (L10); `answer1.txt` (L14); `answer2.txt` (L60) | `maxmean_short_p0109.txt` (L53); `maxmean_middle_p0109.txt` (L98) | `k1; intersection_fit_data` (L49); `k2; intersection_fit_data2` (L94) |
| `sumstat.txt` | `Yamasaki_Go校正後1213_G_log10.txt` (L4) | `SUMSTAT_ct_short_pavlue.txt` (L67); `SUMSTAT_ct_middle_pavlue.txt` (L102) | `answer` (L17); `answer2` (L32); `k1; intersection_fit_data` (L45); `k2; intersection_fit_data2` (L78) |
| `校正前後的繪圖.txt` | No file read. | No table file written. Overwrites plot1 at 12, 39 and 45; contains no explicit image export. | `gene_p_len` (L2) |
| `條狀圖證明.txt` | `short_score.txt` (L9); `mid_score.txt` (L11); `short_mid_score.txt` (L13) | No table file written. Creates three ggplot expressions at 22, 41 and 57; no explicit file export. | No undeclared input object identified in this static pass. |
| `比較分數與無母分析程式.txt` | `all_genes_short_outside_log_noinf.txt` (L5); `all_genes_middle_outside_log_noinf.txt` (L8); `short_mid_key.txt` (L12); `short_mid_all.txt` (L15); `short_CT_G.txt` (L21); `middle_CT_G.txt` (L23) | `mid_result.txt` (L193); `mid_score.txt` (L198); `short_mid_error_bar_ranges.txt` (L202) | `short_point` (L80): Printed before assignment at 95.; `middle_point` (L81): Printed before assignment at 96.; **Fails parse-only check.** |
| `超幾何.txt` | `GO_terms_4899.txt` (L3); `Yamasaki_Go校正後1213_G_log10.txt` (L8); `short_CT_G_score.txt` (L12); `middle_CT_G_score.txt` (L15) | `篩選後大豆ct short0109.txt` (L86) | `answer4` (L113): Inspected but never created.; **Fails parse-only check.** |

## Unresolved upstream connections

- **Yamasaki_Go校正後1213_G_log10.txt**: No active writer for this exact file among these 11 scripts. The correction file and log-transform objects exist as separate steps but are not explicitly connected to this saved file.

- **answer1.txt; answer2.txt**: Hypergeometric code creates answer/answer2 in memory but never actively saves these files. GSEA has one commented pathway write. Folder/version selection is therefore not enforced.

- **CTgenes 含分數 -> short_CT_G_score.txt; middle_CT_G_score.txt; short_CT_G.txt; middle_CT_G.txt**: These scripts read saved CT intermediates; none directly imports the named main CT workbook or defines this extraction/mapping step.

- **WITHOUTINF/outer_data_2.0_short_address.txt; WITHOUTINF/outer_data_2.0_middle_address.txt**: No direct read of either file in these 11 original scripts. Their relationship to corrected score tables or noinf log-score tables requires separate provenance evidence.

- **all_genes_short_outside_log; all_genes_middle_outside_log_noinf.txt; all_genes_short_outside_log_noinf.txt**: Exploratory GSEA relies on an unbound in-memory object; downstream comparisons load saved noinf score files. Their construction is not defined in this folder.

- **SoyNet selected-gene files; crosstalk edge matrices/node lists; plotting summary files**: Loaded as external saved intermediates. No SoyNet/Cytoscape invocation or upstream exporter for these files is present in this folder.

## Per-script evidence and qualifications

### #外部資料匯入與校正繪圖.txt

Map the transformed external expression table to Glyma 2.0; add gene lengths and SNP counts; impute zero SNP counts by length bins; filter mapped rows against GO membership; inspect distributions; apply a correction.

- The external input is the separately saved transformed Pei-Hsiu table, not either WITHOUTINF file (lines 2-4). Its columns 5 and 8 supply p_1 and p_24 (50-51).
- GO-filtered mapped external rows are created at 83-88 and saved at 92. This keeps qualifying rows; it does not deduplicate the external gene identifiers.
- gene_p_len is assigned GO_terms_yamasaki at 104, but p_middle is later assigned out1 at 150; the correction/output branch therefore uses the earlier unfiltered mapped table.
- Only p_middle column 2 is replaced by the power correction at 170. The grid uses Spearman correlation with SNP count at 165.
- The emitted corrected filename lacks the _G_log10 suffix consumed by later statistics; no matching writer for that latter file is present in this folder.

Plot handling: Creates plot1 at 140, but contains no explicit file export.

### crosstalk 20240916.R

Read saved short/middle Jaccard-like edge matrices and node lists; create an interactive combined igraph network.

- Calls install.packages(igraph) at 3; tkplot at 163 requires an interactive Tcl/Tk display.
- All five data files are read explicitly, but their producing code is absent from this script folder.
- Node/edge selection uses manual indices, including loc3=47 at 139 and selected short nodes at 150; these depend on input order and external inspection.

Plot handling: Opens tkplot at 163; contains no explicit saved plot path.

### gsea_菲斯福氣.txt

Exploratory CT-only and full-score GSEA branches with saved CT and pathway intermediates.

- datdir remains the filtered-pathway directory when short_CT_G_score.txt and middle_CT_G_score.txt are read at 36 and 38; both are missing there, but files with these names exist in the CT_genes directory.
- Only middle_esmax.txt is exported (87). The later permutation result objects are not exported.
- The file switches from CTgenes_middle_address for observed scores to all_genes_score for the first permutation branch, and later switches to a pre-existing short score object.
- Saved permutation counts are 100 and 10 (101, 210), although a comment mentions 10000.

Plot handling: None.

### gsea全部.txt

Exploratory short-phase enrichment-score and permutation code using a pre-existing score object.

- answer1.txt is requested from the project root; it is absent at that inferred current location, though alternatives exist under the data directories.
- At 49, read.table(datfile) re-reads answer1.txt; it does not read outfile=short_esmax.txt just constructed at 48.
- The output uses a Linux path /home/hao_yu/soybean/short_esmax.txt (43-45).
- The only pathway-file write is commented out at 122.
- p_t is 2 at 66; this saved exploratory code is not a recorded 10000-permutation run.

Plot handling: None.

### maxmean0127版本.txt

Read corrected scores and saved short/middle pathways; calculate z-scores and the 0127 MaxMean variant and permutations.

- Reads corrected scores plus answer1/answer2, but CT-count columns still depend on prior session objects.
- Writes the 0127 short and middle outputs (70, 133).

Plot handling: None.

### maxmean原版.txt

Read corrected scores and saved short/middle pathways; calculate the older experimental statistic and permutation outputs.

- The three data files are loaded explicitly, but CT-count columns still depend on prior session objects.
- Writes the 0109 short and middle outputs (53, 98).

Plot handling: None.

### sumstat.txt

Read saved corrected scores, sum per-pathway short/middle scores, draw 10000 random gene sets per pathway, BH-adjust and export both phases.

- Reads only Yamasaki_Go校正後1213_G_log10.txt (4); all pathway membership and CT-count objects must exist in the R session.
- Column 2 supplies short scores, column 3 supplies middle scores; it does not load either WITHOUTINF file.
- Simulation comments explicitly say all gene = 42575 (53 and 87), but comments are not proof of the final paper count.

Plot handling: None.

### 校正前後的繪圖.txt

Transform and plot a pre-existing five-column gene table before/after a power correction.

- Reads no file. The expected gene_p_len can be produced by the external import script (104), but no source call declares or enforces this dependency.
- Correction tuning uses Spearman correlation with gene length at 29; the import script instead tunes against SNP count.
- Only column 2 is corrected at 33. Both score columns are transformed to negative base-10 logs at 36.

Plot handling: Overwrites plot1 at 12, 39 and 45; contains no explicit image export.

### 條狀圖證明.txt

Read three manually formatted network-score summary tables and draw grouped bar charts with standard-deviation error bars.

- Calls install.packages for dplyr and ggplot2 (2, 4).
- Reads short_score.txt, mid_score.txt and short_mid_score.txt from the 處理後格式R subdirectory (9, 11, 13).
- Significance-line positions and heights are hard-coded (30-33, 48-51, 64-67), not computed from imported test results.

Plot handling: Creates three ggplot expressions at 22, 41 and 57; no explicit file export.

### 比較分數與無母分析程式.txt

Load separate log-score tables, SoyNet/network-selected genes and CT lists; compare four gene groups using Wilcoxon and t tests and export summary tables.

- Reads six files explicitly (5, 8, 12, 15, 21, 23). Upstream generation of the network gene lists and noinf score tables is not present in these 11 scripts.
- Contains a syntax error at 132, so it cannot be sourced intact.
- Output basenames mid_result.txt and mid_score.txt are used despite combining both phases at 188 and 194; those exact output paths are absent from the current tree.
- The long-format plotting files are in a separate 處理後格式R subdirectory; no automatic reshape/write connection to them is defined here.

Plot handling: None.

### 超幾何.txt

Read GO membership, saved corrected scores and separate short/middle CT lists; select CT-bearing pathways and calculate short-phase hypergeometric results.

- Creates answer, answer2, k1, k2 and intersection indices in memory (42-53); no active write saves answer1.txt or answer2.txt.
- The reference population is unique nonempty all_genes_score column 1 at 63, not the union of CT lists.
- CT lists are flattened to vectors at 13 and 16, then incorrectly indexed as two-dimensional objects at 39 and 50.
- Only a short-phase hypergeometric output is written (86); the middle block is incomplete.
- The file cannot be sourced intact because the un-commented underscore separator at 126 is invalid R syntax.

Plot handling: None.

## Parse-only verification

R 4.5.1 parsed nine scripts successfully. Two fail before execution:

- `超幾何.txt` L126: an un-commented underscore separator is invalid R syntax.
- `比較分數與無母分析程式.txt` L132: an extra comma after a function call is invalid R syntax.

Parsing the remaining files is not an execution test. Hard-coded paths, missing session objects, exploratory fragments, interactive graphics and manual input-format changes remain execution/provenance barriers.

The newly prepared Python/R reanalysis and its completed tests are separate from these original scripts. Successful execution of the new package does not establish that the original research data lineage has been fully reconstructed. No `source()` calls orchestrate the original scripts; upstream saved intermediates and session objects must still be reconciled.
