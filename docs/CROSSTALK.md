# Offline crosstalk before SoyNet

Both languages continue from pathway enrichment through Jaccard crosstalk and gene-list export without contacting SoyNet. The reference workflow uses the author's historical Excel selection: **13 short-term and 11 mid-term pathways**, with two shared GO terms and 22 distinct nodes. It recalculates overlaps from GO annotation rather than copying an old network matrix.

## Run the historical selection

From the repository folder, choose either language:

```sh
python run_pipeline.py --mode full --use-reconstruction-hypotheses --crosstalk reference
```

```sh
Rscript --vanilla run_pipeline.R --mode full --use-reconstruction-hypotheses --crosstalk reference
```

Use `--mode quick` for a 199-draw execution check. Add `--output results/my_run` when the default output exists. The [score-reconstruction caveats](REPRODUCIBILITY.md) still apply.

If enrichment has already completed, reuse its analysis directory without repeating simulations:

```sh
python scripts/run_presoynet.py results/python_full/analysis results/python_crosstalk --selected-pathways data/reference/crosstalk/selected_pathways.tsv --selection-label "Historical Excel selection" --membership annotation --edge-rule historical-code
```

```sh
Rscript --vanilla scripts/run_presoynet.R results/r_full/analysis results/r_crosstalk --selected-pathways data/reference/crosstalk/selected_pathways.tsv --selection-label "Historical Excel selection" --membership annotation --edge-rule historical-code
```

The analysis directory must contain both phases' scores and hypergeometric, SUMSTAT and GSEA tables from the same analysis. Counts, membership sizes, ratios and probability ranges are checked before outputs are written. Input fingerprints and settings accompany each output. Original Excel files and original drive paths are not needed.

## Excel selection and the paper

The author selected the final pathways in Excel using the biological and statistical considerations in the [paper's pathway enrichment section](https://link.springer.com/article/10.1007/s00299-025-03643-2). The paper describes CTgene-count thresholds of 5 for short and 40 for middle, a CTgene ratio of 0.05, and BH significance below 0.005. The reference workflow preserves the actual selection rather than substituting a new list from the reanalysis.

The source is `Table 1. The results of the pathway enrichment analysis with secondary selection..xlsx`, `Short-term!A8:O20` and `Mid-term!B10:J20`. GO IDs were matched by descriptions ignoring case. Both ID sets match the saved February 21 crosstalk lists and local June 2025 manuscript Table S1. The workbook contains values, with no formulas or saved autofilter defining the selection procedure.

The [reference table](../data/reference/crosstalk/selected_pathways.tsv) preserves historical counts and BH values. The [criteria audit](../data/reference/crosstalk/selection_criteria_audit.tsv) evaluates each condition separately:

| Check on selected Excel rows | Short (13) | Middle (11) |
|---|---:|---:|
| Count reaches the phase threshold | 8 | 11 |
| CTgene count / annotation size reaches 0.05 | 6 | 7 |
| Both count and ratio reach their thresholds | 1 | 7 |
| Count or ratio reaches its threshold | 13 | 11 |
| Any main method has BH < 0.005 | 12 | 11 |

An **OR** interpretation of count/ratio is compatible with every selected row. This is an inference from the values, not a recovered formula or proof of the complete selection algorithm. A universal AND rule would lose most selections. Short-term response to abscisic acid has a lowest main-method BH value of 0.009300373, so a universal BH < 0.005 filter would also drop an archived selection. The middle sheet's bold-format footnote mentions 0.0005; the short sheet's copy mentions 0.005. These differences are recorded without changing the historical list.

A separate [programmatic selector](EXCEL_SELECTION.md) now applies explicit rules to the full historical 100/342 candidate tables. Its included candidate configuration reproduces 13/11 exactly, but the original operation sequence remains unconfirmed. The reference workflow continues to preserve the selected Excel IDs.

Three selected short pathways have five annotated genes. They remain in the reference network, with unavailable current enrichment statistics left blank. They are not inserted into the current 6–1,500-gene BH testing family.

## Calculations and outputs

Reference mode uses all unique annotated genes per selected pathway, matching the historical Jaccard calculation. For every unordered pair, `J = intersection size / union size`. Empty padding and self-loops are excluded. All pairs, including zero overlaps, are exported for inspection.

The reference pipeline uses `--edge-rule historical-code`: round the freshly computed coefficient to nine decimals for comparison, then retain values >= 0.004456328 (short) or >= 0.14160401 (middle). These cutoffs come from the recovered original `crosstalk 20240916.R`, which reads the final February 21 matrices. It produces 36 short and 25 middle edges, merged into 61. Both the full-precision coefficient and the filtering value are exported. Rounding is essential for one boundary edge; see the [source trace](../reports/CROSSTALK_TRACE_2026-10-24.md).

The alternative `--edge-rule fraction` ranks positive overlaps and retains `ceil(fraction * positive pairs)`, including cutoff ties. Defaults are 0.8 for short and 0.5 for middle; this interpretation of the paper's fractions produces 67 edges. Standalone commands default to this fraction rule unless explicitly changed. Shared edges are merged with separate phase weights. Isolated selected nodes remain present.

| Output | Contents |
|---|---|
| `*_selection_audit.tsv` | Full current testing family plus extra reference IDs, current statistics and selection provenance |
| `*_selected_pathways.tsv` | Pathways actually used |
| `*_all_pairs.tsv`, `*_jaccard.tsv` | Intersection/union counts, raw and filtering coefficients, retention flags and square matrices |
| `network_nodes.tsv`, `network_edges.tsv` | Merged network suitable for Cytoscape import |
| `pathway_genes.tsv` | Every member with CT and measured-gene flags |
| `network_summary.tsv`, `export_summary.tsv` | Counts and actual edge cutoffs |
| `crosstalk.png` | Deterministic layout grouped by phase; placement is not a clustering result |
| `soynet_inputs/*.txt` | Eight sorted unique gene lists, one ID per line, without headers |

The four export groups are short, middle, shared pathways, and combined. `shared_pathways` contains the union of genes in GO terms selected in both phases, not the intersection of the phase-wide gene sets. Each group has an all-gene file and a CTgene file. Combined CT membership uses both phase CT lists. Isolated pathways contribute genes too. These are handoff groups, not the paper's inferred SoyNet modules.

## Verification and limits

The October 24 check ran both implementations on the saved complete enrichment runs. Nine deterministic tables and all eight gene lists agreed. All 133 within-phase Jaccard pairs matched the saved matrices at printed precision (maximum absolute difference about 5e-10).

| Quantity | Recalculated reference network | Paper |
|---|---:|---:|
| Distinct pathways | 22 | 22 |
| Positive merged edges before filtering | 102 | 102 |
| Edges after the dated original script's rule | 61 | 61 |
| Short CTgenes | 31 | 31 |
| Middle CTgenes | 121 | 122 |
| Combined CTgenes | 145 | 145 |

Pathway identities, overlap coefficients and the final filtered edge count are reproduced using the dated script. The 121/122 middle difference is traced to Glyma.09G149200: the saved middle gene export used 14 pathway sets, including GO:0009793 and GO:0009826, which contain this gene. Neither pathway belongs to the final 11. The 14-set calculation reproduces all 122 archived IDs exactly; the final 11-set calculation remains 121. This is a pathway-version difference, not a lost CT input. Details and the six edges distinguishing the filtering rules are in the [follow-up trace](../reports/CROSSTALK_TRACE_2026-10-24.md).

After generating both reference networks:

```sh
python scripts/verify_presoynet.py results/python_crosstalk results/r_crosstalk --historical-reference
```

Native simulation probabilities are not compared for equality because the languages have different random generators. The comparator checks deterministic downstream outputs, historical overlaps and retained edge identities under the dated script's rule. See the [validation and lineage report](../reports/CROSSTALK_TRACE_2026-10-24.md).

For a separate fresh selection, use `--crosstalk reanalysis`. Defaults are any main-method BH < 0.05, at least one measured CTgene, and measured-gene membership. Standalone commands expose `--alpha`, `--min-ct`, `--min-ratio`, `--method-rule any|all`, `--membership measured|annotation` and edge fractions. This mode is labelled `bh_reanalysis` and does not claim to recover the Excel choices.

SoyNet queries, inferred gene networks, topology-based gene prioritization and experimental validation remain outside this offline step.
