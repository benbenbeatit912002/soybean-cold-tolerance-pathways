"""Independently compare fresh outputs with original files, never alter outputs.

Usage: python scripts/verify_short_recalculation.py ORIGINAL_PROJECT RESULTS_ROOT
Expects python/, r/, python_hypothesis/, and r_hypothesis/ result directories.
The report separates numerical agreement from historical provenance.
"""
import argparse
from collections import defaultdict, deque
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"python"))
from ctpathway import canonical_id, read_tsv, write_tsv
from historical_short import SOURCE_FILES
from prepare_historical_correction import REFERENCE_HASHES


def compare(source_root, results_root):
    source = source_root/"研究法/大豆基因/數據/外部基因校正"
    hashes = {name: hashlib.sha256((source/name).read_bytes()).hexdigest() for name in REFERENCE_HASHES}
    if hashes != REFERENCE_HASHES:
        raise ValueError("Original source changed: re-audit before comparing")
    inputs = read_tsv(source/"gene_p_len.txt")
    corrected = read_tsv(source/"Yamasaki_Go校正後1213.txt")
    logged = read_tsv(source/"Yamasaki_Go校正後1213_G_log10.txt")
    groups = defaultdict(deque)
    for row in logged:
        groups[(canonical_id(row["gene_2.0"]), int(row["len"]))].append(row)
    expected = []
    for raw, adjusted in zip(inputs, corrected):
        key = (canonical_id(raw["gene_2.0"]), int(raw["len"]))
        saved = groups[key].popleft()
        if key != (canonical_id(adjusted["gene_2.0"]), int(adjusted["len"])):
            raise ValueError("Saved raw/corrected order changed")
        if abs(-math.log10(float(adjusted["p_1"]))-float(saved["p_1"])) > 5e-8:
            raise ValueError("Invalid occurrence join")
        expected.append(saved)
    if any(groups.values()) or len(inputs) != len(corrected):
        raise ValueError("Source record count mismatch")

    report = {"source_sha256": hashes,
              "upstream_source_sha256": {name: hashlib.sha256((source_root/path).read_bytes()).hexdigest()
                                         for name, path in SOURCE_FILES.items()},
              "records": len(inputs), "runs": {},
              "scope": "Recalculation from four upstream files; saved scores used here only for verification",
              "historical_zero_repair_provenance_confirmed": False,
              "all_paper_methods_and_pathway_results_reproduced": False}
    outputs = {}
    for mode in ("python", "r", "python_hypothesis", "r_hypothesis"):
        rows = read_tsv(results_root/mode/"recalculated_short.tsv")
        outputs[mode] = rows
        if len(rows) != len(inputs):
            raise ValueError(f"Wrong record count: {mode}")
        differences = []
        input_mismatches = []
        for index, (row, raw, saved) in enumerate(zip(rows, inputs, expected), start=2):
            identity_matches = (int(row["source_row"]) == index and
                                row["gene_id"] == canonical_id(raw["gene_2.0"]))
            input_matches = all(float(row[key]) == float(raw[old]) for key, old in
                                [("raw_short", "p_1"), ("raw_middle", "p_24"),
                                 ("gene_length", "len"), ("historical_snp", "snp")])
            if not identity_matches or not input_matches:
                input_mismatches.append(index)
            actual = Decimal(row["score_short"])
            original = Decimal(saved["p_1"])
            if not actual.is_finite() or actual != original:
                differences.append({"source_row": index, "gene_id": row["gene_id"],
                                    "calculated_score": str(actual), "saved_score": str(original)})
        write_tsv(differences, results_root/mode/"verification_mismatches.tsv",
                  ["source_row", "gene_id", "calculated_score", "saved_score"])
        grids = read_tsv(results_root/mode/"short_lambda_grid.tsv")
        best = min(grids, key=lambda row: float(row["absolute_spearman"]))
        report["runs"][mode] = {"records": len(rows), "input_records_matching_exactly": len(rows)-len(input_mismatches),
            "output_sha256": hashlib.sha256((results_root/mode/"recalculated_short.tsv").read_bytes()).hexdigest(),
            "input_mismatch_source_rows": input_mismatches, "short_lambda": float(best["lambda"]),
            "score_values_matching_exactly": len(rows)-len(differences), "score_mismatches": differences,
            "hypothesis_enabled": mode.endswith("hypothesis"),
            "comparison": "exact decimal equality at saved output precision, without a tolerance"}
        if mode.endswith("hypothesis"):
            report["runs"][mode]["repair_audit"] = read_tsv(results_root/mode/"zero_repair_hypothesis.tsv")
    for suffix in ("", "_hypothesis"):
        left, right = outputs["python"+suffix], outputs["r"+suffix]
        report["cross_language"+suffix] = {
            "exported_score_values_equal": sum(Decimal(a["score_short"]) == Decimal(b["score_short"]) for a,b in zip(left,right)),
            "maximum_adjusted_probability_difference": max(abs(float(a["adjusted_short"])-float(b["adjusted_short"])) for a,b in zip(left,right))}
    report["numerical_hypothesis_reproduces_all_saved_short_scores"] = all(
        not report["runs"][mode]["score_mismatches"] and not report["runs"][mode]["input_mismatch_source_rows"]
        for mode in ("python_hypothesis", "r_hypothesis"))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("results_root", type=Path)
    args = parser.parse_args()
    report = compare(args.source_root, args.results_root)
    destination = ROOT/"reports/short_recalculation_verification_2026-10-21.json"
    destination.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if not report["numerical_hypothesis_reproduces_all_saved_short_scores"]:
        sys.exit(1)
