# Data, provenance, and reuse

The recommended Python and R pipelines use only included files. Upstream source files are preserved byte for byte. Original filenames are retained in manifests as provenance, not required directory names.

## Files used by the main workflow

| Packaged path | Source role |
|---|---|
| `upstream/expression.tsv` | Converted Peixiu expression table; original columns 3, 5, and 8 provide source gene ID, short p-value, and middle p-value |
| `input/gene_correspondence.tsv` | Original correspondence table; first two columns are aliases, third is the target gene ID |
| `upstream/gene_lengths.tsv` | Original gene-length table; columns 1 and 5 give gene ID and length in base pairs |
| `upstream/snp_markers.txt` | Original gene blocks with AX markers and END delimiters; empty blocks preserve zero counts |
| `input/ct_short.tsv`, `input/ct_middle.tsv` | Primary CTgenes membership and prior scores |
| `input/pathways.gmt` | Prepared GO memberships from GO_terms_4899.txt |

`upstream/inputs.tsv` lists all four reconstruction sources, paths relative to its own directory, original relative names, byte counts, SHA-256, and MD5. The correspondence file is shared with `input/` to avoid a second copy. Python verifies SHA-256 and MD5; base R verifies MD5. A changed or missing packaged source is an error, never a reason to search another folder.

`.gitattributes` disables Git line-ending conversion for data and archived scripts to preserve hashes. Do not open and resave source files in Excel. Create a separate processed table if you need another transformation.

## Other included data

- `input/external_short.tsv` and `external_middle.tsv`: prepared WITHOUTINF data for the older alternative analyses. Column 5 is fold change; column 6's transformation remains unresolved. These are not inputs to the new main score-recalculation workflow.
- `reference/`: historical inputs and expected short scores for independent checks. The calculation and simulations never use expected scores.
- `../examples/data/`: invented demonstration data with TOY gene IDs, not research findings.

The source manifests in `../reports/` record the original staging of CT lists, GO memberships, external files, and archival scripts. Detailed fields and exclusions are in the [data dictionary](../docs/DATA_DICTIONARY.md).

## Source boundary

The expression input comes from the supplied converted TXT table. Yamasaki Excel files are not runtime inputs, and their full conversion chain to this TXT remains unresolved. See [reproducibility status](../docs/REPRODUCIBILITY.md).

## Reuse

The analysis accompanies Liu et al. (2025), DOI 10.1007/s00299-025-03643-2. Collaborator-provided material is attributed in the provenance documentation. A software license and data redistribution license have not yet been assigned. Including inputs here does not establish their original licensing terms.

Software and data licenses remain unassigned. Data provenance is documented above; provenance does not establish redistribution or reuse permissions. Historical Chinese filenames and comments remain unchanged in archival material.
