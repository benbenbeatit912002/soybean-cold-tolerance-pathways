"""Stage the verified historical correction input without changing its records.

Usage: python scripts/prepare_historical_correction.py SOURCE_ROOT OUTPUT.tsv
This file already contains earlier SNP imputations. Only its remaining zero is
handled by the correction CLI; this is not a reconstruction from raw SNP markers.
"""
import argparse
from collections import defaultdict, deque
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
from ctpathway import canonical_id, read_tsv, write_tsv


# These are independently audited expected-result snapshots for verification.
# They must not be read by a score calculator. A changed source needs a fresh audit.
REFERENCE_HASHES = {
    "gene_p_len.txt": "3ff5d73fc5c42a542a0eb7765b733086a2022eb21b33ff6f84b40e42da467be9",
    "Yamasaki_Go校正後1213.txt": "1d6818f1668afc1875d9d8570cfb6c4c479a0207ed942ba2951dbfb69d485e29",
    "Yamasaki_Go校正後1213_G_log10.txt": "f476cf50e2212bbe0e4f6d6bcd0addd2809db04e7c77ab6d4b2588217fa4f9e5"}


def stage_short_reference(source_directory, rows, destination):
    """Preserve saved decimals and identify duplicate records by occurrence."""
    for name, expected in REFERENCE_HASHES.items():
        if hashlib.sha256((source_directory/name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Historical reference source changed; audit it before restaging: {name}")
    corrected = read_tsv(source_directory/"Yamasaki_Go校正後1213.txt")
    logged = read_tsv(source_directory/"Yamasaki_Go校正後1213_G_log10.txt")
    if not len(rows) == len(corrected) == len(logged):
        raise ValueError("Historical source record counts differ")
    groups = defaultdict(deque)
    for line, row in enumerate(logged, start=2):
        groups[(canonical_id(row["gene_2.0"]), int(row["len"]))].append((line, row))
    reference = []
    largest_log_difference = 0.0
    for row, old in zip(rows, corrected):
        key = (canonical_id(row["gene_id"]), int(row["gene_length"]))
        if key != (canonical_id(old["gene_2.0"]), int(old["len"])) or not groups[key]:
            raise ValueError("Historical record identity differs")
        line, saved = groups[key].popleft()
        difference = abs(-math.log10(float(old["p_1"]))-float(saved["p_1"]))
        if difference > 5e-8 or float(old["snp"]) != float(saved["snp"]):
            raise ValueError("Occurrence join failed independent score/SNP validation")
        largest_log_difference = max(largest_log_difference, difference)
        reference.append({key: row[key] for key in
                          ("source_row", "gene_id", "raw_short", "gene_length", "snp")})
        reference[-1].update(reference_adjusted_short=old["p_1"],
                             reference_score_short=saved["p_1"],
                             reference_snp_short=saved["snp"], reference_log_row=line)
    if any(groups.values()):
        raise ValueError("Unmatched historical logged records")
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(reference[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(reference)
    content = buffer.getvalue().encode("utf-8")
    if destination.exists() and destination.read_bytes() != content:
        raise ValueError("Refusing to replace a different saved reference; choose a new destination")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(content)
    provenance = {
        "author_designation": "Original short-term scores designated as paper results on 2026-10-21",
        "source_sha256": REFERENCE_HASHES, "reference_sha256": hashlib.sha256(content).hexdigest(),
        "records": len(reference), "unique_genes": len({canonical_id(r["gene_id"]) for r in rows}),
        "join": "canonical gene ID, gene length, and occurrence order within each pair",
        "short_score_source": "p_1 text in Yamasaki_Go校正後1213_G_log10.txt",
        "adjusted_short_source": "p_1 text in Yamasaki_Go校正後1213.txt",
        "independent_log_check_max_absolute_difference": largest_log_difference,
        "policy": "Preserve saved strings exactly; no averaging, recalculation, or inferred exception formula."}
    destination.with_suffix(".provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(f"Preserved {len(reference)} original short scores: {destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--short-reference", type=Path,
                        help="also freeze the audited original short scores at this path")
    args = parser.parse_args()
    source = args.source_root/"研究法/大豆基因/數據/外部基因校正/gene_p_len.txt"
    source_rows = read_tsv(source)
    rows = [{"source_row": i+2, "gene_id": row["gene_2.0"], "raw_short": row["p_1"],
             "raw_middle": row["p_24"], "gene_length": row["len"], "snp": row["snp"]}
            for i, row in enumerate(source_rows)]
    if args.short_reference:
        stage_short_reference(source.parent, rows, args.short_reference)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(rows, args.output)
    metadata = {"source": str(source.resolve()),
                "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "rows": len(rows), "zero_snp_rows": sum(float(r["snp"]) == 0 for r in rows),
                "policy": "Preserve all rows, scores, lengths, and stored SNP counts; rename columns only.",
                "scope": "Saved input already contains historical imputations; not raw SNP reconstruction."}
    args.output.with_suffix(".provenance.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(f"Staged {len(rows)} historical records: {args.output}")


if __name__ == "__main__":
    main()
