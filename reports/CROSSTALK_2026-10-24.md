# Offline crosstalk verification — 2026-10-24

**Follow-up:** the [dated-script and Excel trace](CROSSTALK_TRACE_2026-10-24.md) resolves the edge-count procedure and identifies the middle CTgene difference. This report retains the initial fraction-rule results as a historical checkpoint.

The extension implements the workflow up to SoyNet in independent Python and base-R versions. Existing enrichment calculations are unchanged. The author's manually selected Excel pathways are preserved as historical references, while overlaps and gene exports are recalculated from packaged annotation and CT lists.

## Source checks

- `Short-term!A8:O20` contains 13 pathways; `Mid-term!B10:J20` contains 11 in the secondary-selection workbook.
- Both ID sets match the saved February 21 crosstalk lists and the local June 2025 manuscript Table S1.
- The workbook has no formula cells or saved autofilter. The paper supplies the selection considerations; individual Excel rows are assessed in `data/reference/crosstalk/selection_criteria_audit.tsv`.
- All selected rows meet count OR ratio thresholds; count AND ratio retains only 1 short and 7 middle rows. OR is compatible with the selected values, not a recovered formula or proof that it selects only these rows from the original candidate population.
- One selected short row, GO:0009737, has no main-method BH value below 0.005. Its lowest value is 0.009300373.
- Exact source filenames, ranges and SHA-256 checksums are in `data/reference/crosstalk/provenance.json`. Identity with the final publisher supplement has not been verified.

## Independent execution

Both implementations were run against the complete Python/R enrichment outputs saved on October 21, using all annotated pathway genes and the same historical reference list. Existing 10,000-draw enrichment results were reused; this downstream change does not require new score simulations.

The cross-language comparison passed for nine deterministic tables and all eight exported gene lists. All 133 unordered within-phase Jaccard pairs matched the independently saved matrices within 5.1e-10; the largest observed difference was 5.000000413701855e-10.

| Phase | Selected pathways | Positive pairs | Retained pairs | Cutoff |
|---|---:|---:|---:|---:|
| Short | 13 | 48 | 39 | 0.0040032025620496394 |
| Middle | 11 | 55 | 28 | 0.12783906554185595 |

The union contains 22 nodes and 102 positive edges before filtering. The explicit top-fraction rule retains 67 merged edges. Retention ranks only positive undirected pairs, rounds the requested count upward and includes cutoff ties. The paper's 61-edge filtered graph is not reproduced by this operational definition. Earlier scripts use manual cutoffs and a different middle-pathway list; they do not supply a confirmed complete specification of the final filtering procedure.

| Export group | Unique genes | CTgenes |
|---|---:|---:|
| Short | 4,697 | 31 |
| Middle | 3,869 | 121 |
| Shared pathways | 1,863 | 68 |
| Combined | 6,102 | 145 |

The middle CTgene total differs from the paper's 122. The current membership is recalculated, not overwritten with the Excel counts. The other displayed CTgene totals agree with the corresponding paper totals. The exported groups do not claim to reconstruct the five SoyNet modules.

## Validation scope

Known-set tests cover exact intersections/unions, zero overlaps, self-loop exclusion, ties, isolates, empty and single-node graphs, fraction boundaries, explicit reference selection and rejection of mixed metadata or invalid probabilities. Portable entry-point tests execute both languages from another directory with spaces and check refusal to overwrite existing results. Package checks additionally rerun both real-data quick pipelines with reference crosstalk and compare their outputs against the archived coefficients.

The statistical core is unchanged. Native simulated probabilities may differ between languages; their equality is not required for the historical-list network. New reanalysis selection remains a separate labelled option. No SoyNet service was contacted and no experimental validation was run.

Commands and output definitions: [Crosstalk guide](../docs/CROSSTALK.md).
