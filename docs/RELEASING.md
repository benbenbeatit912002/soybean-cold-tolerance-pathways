# Prepare and validate a release

Build a source package, validate it from a clean extraction, and retain its checksum with the validation report.

## Build a source package

From the repository root:

```sh
python scripts/package_project.py
```

The output is `dist/soybean-cold-tolerance-pathways-source.zip` with a companion SHA-256 file. It includes code, all required input data, examples, documentation, and curated reports. It excludes Git metadata, virtual environments, caches, test scratch files, generated results, and previous package-test extractions.

This source package is the download to test. It must not contain precomputed results that could hide missing steps.

## Test the extracted package

Install the Python requirements and R first. Then:

```sh
python scripts/check_package.py
```

If Rscript is not on PATH:

```powershell
python scripts/check_package.py --rscript "C:/Program Files/R/R-4.5.1/bin/Rscript.exe"
```

The checker extracts to a new directory containing spaces, runs Python and R tests, checks shared-null agreement, runs both demos and real-data quick workflows, checks that strict unresolved-score runs stop, and independently verifies the recalculated outputs. The entry-point commands start outside the extracted repository. No original project folder or saved R session is needed.

For a complete additional 10,000-draw rerun in both languages, including experimental MaxMean:

```sh
python scripts/check_package.py --full
```

Add `--rscript` as above when needed. Expect this to take longer than the default check. Logs and extracted outputs remain under `results/package_checks/`. The report is `reports/latest_portable_package_check.json`; it records the exact tested ZIP hash, file hashes, exit codes, durations, and scope. A report is evidence only for that tested snapshot. Later documentation-only packaging updates should be identified separately.

## Optional results bundle

```sh
python scripts/package_project.py --include-results
```

This creates a separate `*-with-results.zip`. It excludes nested package-test extractions to avoid recursive bundles. Curate the results you intend to distribute before publishing; generated data are not committed by default.

## Publish a reviewed version

Use Git to track source changes and a version tag to identify a reviewed snapshot. A first documented reanalysis can be labelled `v0.1.0`; the tag is not created by these scripts. Attach the source package, checksum, validation summary, and any curated results archive to a GitHub Release.

Before publication, select the code license and record data reuse terms. Preserve `CITATION.cff` and the explicit reconstruction caveats. Do not describe unresolved pathway differences as exact paper reproduction. Optional Zenodo integration can archive releases for research citation.

The CI configuration runs statistical, demo, and real-data quick checks on Windows and Linux. Local validation does not imply those remote jobs have run.
