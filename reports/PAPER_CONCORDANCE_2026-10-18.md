# Paper concordance audit — 18 October 2026

**Verdict: execution is verified; reproduction of the published analysis is not verified.** The current Python and R workflows are a reanalysis of the WITHOUTINF inputs. The audit identifies substantive differences in score preparation, the gene universe, pathway eligibility, and final selection. Agreement between the two implementations cannot validate their shared scientific assumptions.

**Follow-up evidence:** the [historical data-flow audit](HISTORICAL_DATA_FLOW_2026-10-18.md) now reconstructs the original ID mapping, identifies an intermediate 42,576-record count before one empty-target exclusion, and numerically reproduces both saved FOSCO phases for all 42,574 positive-SNP rows. It also confirms cached `10^logFC` columns and explains most logged-score discrepancies through rounding. These resolve several earlier provenance uncertainties; the new WITHOUTINF workflow and final paper results remain separate reproduction questions.

This audit examines repository commit `5691efb265afcc8e357d56272964abd39691f4a7`, the saved full-run outputs, the published article, and local historical files. It does not change numerical code, configurations, inputs, or results.

## Evidence and limits

The [published article](https://link.springer.com/article/10.1007/s00299-025-03643-2) reports 170 CT genes, 4,896 GO terms, and a validation universe of 42,576 genes. It describes a reciprocal-symmetric fold-change transformation followed by FOSCO. Its stated pathway bounds are 6–1,500 genes with CT overlap, yielding 100 short and 342 middle candidates. It reports 13 short and 11 middle selected pathways, with 31 and 122 CT genes, respectively. SUMSTAT uses 10,000 random gene sets. The prose also describes secondary count, ratio, and significance criteria.

The article text was accessible. **The final publisher-hosted supplementary workbooks were not successfully retrieved:** download attempts returned a client-challenge HTML page or a network error. Local manuscript workbooks and the June 2025 submission supplement were read, but their identity with the final published supplements is unverified. Differences below must not be attributed to the final published tables without that check.

Local paths in this report are relative to the author's original `D:/Bioinformatics project` folder unless they identify files within this repository. The original workbooks were opened read-only. Source workbook coordinates are included for review.

## What has actually passed

The [full simulation validation](FULL_SIMULATION_VALIDATION.md) documents completed Python and R runs for columns 5 and 6, both phases, and SUMSTAT, positive-tail GSEA, and experimental MaxMean. Each uses 10,000 random gene sets per distinct measured pathway size. Exact hypergeometric probabilities, BH adjustment, deterministic statistics, and shared-null fixtures were checked. Real-data FOSCO was not run from a confirmed historical input; its numerical implementation and command-line workflow were tested with synthetic data.

These checks establish successful execution and consistency under the recorded configurations. They do not establish that the input score, background, or selection rules match the published analysis. The raw-column-5 runs are especially unsuitable as a paper-reproduction claim in light of the provenance recovered below.

## 1. Input identity and score meaning

The primary lists contain 44 short and 143 middle genes, with 170 unique genes in their union. Their IDs agree exactly with the local June `TableS4_CTgenes candidates (Kao et al., 2022).xlsx`; short scores match exactly and middle scores differ only by rounding, at most `4.98e-9`.

`GO_terms_4899.txt` contains 4,896 terms and 44,964 unique annotated genes. Its declared sizes match its membership lists. It contains GO membership, not numerical gene scores.

Column meanings can now be traced beyond the author's description of scores supplied by Edward Kao:

| External column | Recovered evidence | Consequence |
|---|---|---|
| 2–4 | Historical headers identify three treatment replicates. | These are expression measurements, not pathway statistics. |
| 5 | Historical headers identify 1 hr / 24 hr fold change. The original workbook confirms treatment-replicate total divided by control-replicate total. | Direct larger-is-stronger use does not implement reciprocal-symmetric change magnitude. |
| 6 | Historical header says `Pvalue`, but values exceed one and the workbook contains saved values without formulas. | Its mathematical construction and historical role remain unresolved; it cannot be treated as a raw probability. |

Header sources, line 1:

- `研究法/大豆基因/數據/外部資料匯入/舊版/耐寒資料庫 short _pvalue.txt`
- `研究法/大豆基因/數據/外部資料匯入/舊版/耐寒資料庫 middle _pvalue.txt`

Workbook source: `研究法/大豆基因/數據/外部資料匯入/Pvalue轉換_Yamasaki (2016) 耐寒資料庫0 1 24hr.xlsx`, sheets `1hr cold treatment` and `24hr cold treatment`, headers `A1:F1`. All 40,061 short and 40,332 middle WITHOUTINF five-number tuples match the corresponding treatment worksheet within absolute/relative tolerance `1e-8`. In each original worksheet, all 45,173 finite fold changes with nonzero control totals agree with the treatment/control calculation within `1e-12`, matching by original gene ID. Both treatment worksheets contain zero formula elements in their XML.

The current full column-5 configurations set `score_kind="score"`, so they leave fold change unchanged. For example, FC=2 and FC=0.5 have identical absolute log magnitude but different raw scores. This changes SUMSTAT totals and GSEA rankings. Selecting the existing `fold_change` option would address the magnitude transformation only; it would not reconstruct real-data FOSCO or establish the historical duplicate policy.

Column 6 is neither `abs(ln(FC))` nor the negative base-10 log of the separate saved raw-p branch in the comparisons performed. That branch is `耐寒資料庫0 1 24hr(轉換後)培修版.txt`. The archived correction script loads it at `legacy/#外部資料匯入與校正繪圖.txt:3–4` and extracts columns 5 and 8 at lines 50–51. Consequently, a shared informal label such as “score” or “Pvalue” is insufficient to identify the intended correction input.

## 2. Background and retained CT membership differ

| Quantity | Short | Middle |
|---|---:|---:|
| Primary CT genes | 44 | 143 |
| Raw WITHOUTINF rows | 40,061 | 40,332 |
| Unique resolved external genes in the reanalysis | 39,958 | 40,235 |
| CT genes retained in that measured universe | 29 | 138 |
| Current tested pathways | 2,638 | 2,637 |
| Current tested pathways with CT overlap | 96 | 335 |

The legacy statistics use another branch: `研究法/大豆基因/數據/外部基因校正/Yamasaki_Go校正後1213_G_log10.txt`. It has 42,575 rows but 42,567 unique genes, including eight conflicting duplicated IDs. The legacy short CT score file contains 41 genes. See `legacy/sumstat.txt:2–4`, `legacy/超幾何.txt:6–15`, and [input_audit.json](input_audit.json).

Neither the current resolved universes nor that saved historical file directly establish the reported 42,576-gene universe. This matters for the hypergeometric denominator, the CT foreground, pathway sizes, and the random-gene null distributions. Mean aggregation of duplicate scores is an explicit reanalysis assumption, not a recovered historical decision.

## 3. Pathway size bounds and the multiple-testing family differ

Independent recomputation from the staged full annotation and primary CT lists gives:

| Eligibility rule, requiring at least one primary CT gene | Short | Middle |
|---|---:|---:|
| Full annotation size 6 through 1,500 | 97 | 337 |
| Historical code bounds 5 through 1,499 | 100 | 342 |

The entire difference is three short and five middle size-5 pathways. There are no size-1,500 cases offsetting the difference. The historical bounds reproduce the reported candidate counts, whereas the stated bounds do not reproduce those counts on the supplied annotation. This is evidence of a reporting/version discrepancy to resolve, not permission to silently select whichever setting matches a count.

The current implementation instead filters on size after intersection with each measured universe and applies BH across all eligible pathways, including those without CT overlap. That is a different family from CT-containing candidates selected on full annotation size. Changing the family can change adjusted p-values even if a pathway's raw statistic stays unchanged.

The three short size-5 terms are `GO:0009415`, `GO:0010371`, and `GO:0035264`. All occur in the local selected-pathway table and are excluded by the current minimum measured size of 6.

## 4. The local selected pathways are not fully reproduced

Source: `研究法/大豆基因/期刊撰寫/table/Table 1. The results of the pathway enrichment analysis with secondary selection..xlsx`.

- `Short-term!A8:C20`: 13 selected pathways.
- `Mid-term!B10:D20`: 11 selected pathways.
- The combined list contains 22 distinct GO IDs.
- Using full annotation membership and the primary CT lists, their CT unions contain 31 short and 121 middle genes. The middle count does not establish the reported 122.
- Restricting those same selected memberships to the current external universes leaves 20 short and 119 middle CT genes. These are membership checks of the local selected list, not a newly selected significant result.
- The local short table lists 14 CT genes for `GO:0045893` and 11 for `GO:0050832`; the supplied full annotation and primary CT lists give 15 and 12. The other selected short counts and all selected middle counts agree.

Per-pathway coordinates, membership counts, and the source workbook SHA-256 are saved in [paper_concordance_counts_2026-10-18.json](paper_concordance_counts_2026-10-18.json).

## 5. Final selection cannot yet be reconstructed unambiguously

The local June `TableS1.docx` has Hypergeometric, SUMSTAT, and GSEA columns. It does not define an executable common AND/OR rule for CT counts, CT ratios, nominal p-values, BH values, and a top-fraction cutoff. Its short footnote uses nominal `<0.005`; its middle footnote uses `<0.0005`. Neither defines a separate ranking denominator or rounding rule for the top fraction.

The selected rows also prevent assuming that every stated cutoff applies jointly to every method:

- Short response to water has size 5 and one CT gene; two other short rows also have size 5.
- Short response to cold has a CT ratio of `12/1113`, about 0.0108.
- Middle protein targeting to membrane has a CT ratio of `40/1016`, about 0.0394.
- In the short table, response to water has bold raw hypergeometric p=0.00481 and bold BH q=0.0283. Photoprotection has bold SUMSTAT p=0.0003 and q=0.0075. Formatting does not establish the unspecified BH threshold.

DOCX OOXML coordinates: short Table 1, rows 5, 9, 11, 13 and footnote row 16; middle Table 2, row 3 and footnote row 14. These are local submission evidence, not verified final-publisher evidence. The current package does not implement the paper's final secondary-selection stage.

## 6. Local result versions disagree

The June `TableS2_Pathway enrichment analysis.xlsx` and the local Table 1 workbook contain different short-term hypergeometric results:

| GO ID | June TableS2 location, sheet `Hypergeometric(Short-term)` | TableS2 raw p | Local Table 1 raw p |
|---|---|---:|---:|
| GO:0009409 | B35:G35 | 0.01219794 | 3.78e-10 |
| GO:0009631 | B45:G45 | 2.18e-7 | 8.39e-12 |
| GO:0045893 | B100:G100 | 0.000400304 | 7.57e-13 |
| GO:0050832 | B109:G109 | 0.003002138 | 8.77e-10 |

For the last two rows the CT counts also differ, as described above. The local Table 1 cold-response and cold-acclimation probabilities agree with a hypergeometric universe of 42,567 and foreground of 41. This observation does not establish those parameters for TableS2 or for the final published analysis.

June supplement folder: `研究法/大豆基因/期刊撰寫/期刊要求的檔案/20250623 plant cell report/Supplementary data`.

## Required work before claiming paper reproduction

1. Identify the final published supplementary versions and their actual source inputs; reconcile the conflicting local tables.
2. Reconstruct the intended external score and real-data FOSCO chain, including SNP mapping, missing values, lambda search, and duplicate handling. Recover column 6's formula if that column is used.
3. Fix a documented gene universe and CT foreground for each method, with an explicit gene-by-gene explanation for exclusions.
4. Resolve size bounds, annotation versus measured membership, the BH family, GSEA direction/weight, and the method-specific secondary-selection rule.
5. Run Python and R using those confirmed settings and compare GO IDs, observed statistics, raw/adjusted p-values, and selected gene membership with the final supplementary tables. Monte Carlo comparisons need sampling uncertainty, not exact equality across native random generators.

Until those steps are complete, the accurate description is **a runnable, tested reanalysis with unresolved differences from the published workflow**. No new arithmetic defect was established by this audit; scientific concordance is nevertheless not passed.
