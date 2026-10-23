"""Check that both implementations retain every original saved short score.

Usage: python scripts/verify_paper_short_reference.py SOURCE_ROOT RUN_DIRECTORY
Recomputed FOSCO values are checked separately, never substituted for references.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"python"))
import numpy as np
from ctpathway import canonical_id, read_tsv, write_tsv


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_fosco(p, count, lam):
    return 1.0 if p == 1 else -math.expm1(count**lam * math.log1p(-p))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--input", type=Path, default=ROOT/"data/reference/historical_input.tsv")
    parser.add_argument("--reference", type=Path, default=ROOT/"data/reference/paper_short.tsv")
    args = parser.parse_args()
    run = args.run_directory
    input_rows, reference = read_tsv(args.input), read_tsv(args.reference)
    provenance = json.loads(args.reference.with_suffix(".provenance.json").read_text(encoding="utf-8"))
    directory = args.source_root/"研究法/大豆基因/數據/外部基因校正"
    for name, expected in provenance["source_sha256"].items():
        assert sha256(directory/name) == expected, f"Changed source: {name}"
    assert sha256(args.reference) == provenance["reference_sha256"]
    original = read_tsv(directory/"gene_p_len.txt")
    old_adjusted = read_tsv(directory/"Yamasaki_Go校正後1213.txt")
    old_logged = read_tsv(directory/"Yamasaki_Go校正後1213_G_log10.txt")
    outputs = {lang: read_tsv(run/lang/"corrected_scores.tsv") for lang in ("python", "r")}
    assert len(input_rows) == len(reference) == len(original) == len(old_adjusted) == len(old_logged)
    assert all(len(rows) == len(original) for rows in outputs.values())

    # Physical row numbers disambiguate repeated genes. Verify the reference
    # against original files, not only against the staged reference itself.
    source_rows = list(range(2, len(original)+2))
    assert [int(row["source_row"]) for row in input_rows] == source_rows
    assert [int(row["source_row"]) for row in reference] == source_rows
    assert sorted(int(row["reference_log_row"]) for row in reference) == source_rows
    for source, staged, saved, corrected in zip(original, input_rows, reference, old_adjusted):
        logged = old_logged[int(saved["reference_log_row"])-2]
        gene = canonical_id(source["gene_2.0"])
        assert gene == canonical_id(staged["gene_id"]) == canonical_id(saved["gene_id"])
        assert gene == canonical_id(corrected["gene_2.0"]) == canonical_id(logged["gene_2.0"])
        for field, old_field in (("raw_short", "p_1"), ("gene_length", "len"), ("snp", "snp")):
            assert float(staged[field]) == float(saved[field]) == float(source[old_field])
        assert float(staged["raw_middle"]) == float(source["p_24"])
        assert int(source["len"]) == int(corrected["len"]) == int(logged["len"])
        assert saved["reference_score_short"] == logged["p_1"]
        assert saved["reference_adjusted_short"] == corrected["p_1"]
        assert float(saved["reference_snp_short"]) == float(logged["snp"]) == float(corrected["snp"])
        assert abs(-math.log10(float(corrected["p_1"]))-float(logged["p_1"])) <= 5e-8

    counts = [float(row["snp"]) for row in input_rows]
    bins = [(int(row["gene_length"])-1)//2000 for row in input_rows]
    donors = [n for n in counts if n > 0]
    expected_counts = counts.copy()
    for i, count in enumerate(counts):
        if count == 0:
            group = [n for n, b in zip(counts, bins) if n > 0 and b == bins[i]] or donors
            expected_counts[i] = round(math.fsum(group)/len(group))

    summaries = {lang: {row["phase"]: row for row in read_tsv(run/lang/"correction_summary.tsv")}
                 for lang in outputs}
    for lang, rows in outputs.items():
        assert [int(row["source_row"]) for row in rows] == source_rows
        assert [float(row["snp"]) for row in rows] == expected_counts
        assert [float(row["snp_original"]) for row in rows] == counts
        for row, saved, staged in zip(rows, reference, input_rows):
            assert row["gene_id"] == canonical_id(staged["gene_id"])
            assert int(row["gene_length"]) == int(staged["gene_length"])
            assert row["score_short"] == saved["reference_score_short"], "Saved score text changed"
            assert row["adjusted_short"] == saved["reference_adjusted_short"]
            assert row["short_score_source"] == "historical_reference"
            assert int(row["reference_log_row"]) == int(saved["reference_log_row"])
            assert float(row["reference_snp_short"]) == float(saved["reference_snp_short"])
            for phase in ("short", "middle"):
                assert float(row[f"raw_{phase}"]) == float(staged[f"raw_{phase}"])
                lam = float(summaries[lang][phase]["lambda"])
                expected = stable_fosco(float(staged[f"raw_{phase}"]), float(row["snp"]), lam)
                prefix = "recomputed_" if phase == "short" else ""
                assert math.isclose(float(row[f"{prefix}adjusted_{phase}"]), expected, rel_tol=1e-12, abs_tol=0)
                assert math.isclose(float(row[f"{prefix}score_{phase}"]), -math.log10(expected), rel_tol=1e-12, abs_tol=1e-12)
        assert summaries[lang]["short"]["score_source"] == "historical_reference"
        assert summaries[lang]["middle"]["score_source"] == "recomputed_fosco"

    language_errors = {}
    for column in ("score_short", "adjusted_short", "recomputed_score_short",
                   "recomputed_adjusted_short", "score_middle", "adjusted_middle"):
        a, b = [np.array([float(row[column]) for row in outputs[lang]]) for lang in outputs]
        assert np.isfinite(a).all() and np.isfinite(b).all()
        np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12)
        language_errors[column] = float(np.max(np.abs(a-b)))
    for phase in ("short", "middle"):
        assert summaries["python"][phase]["lambda"] == summaries["r"][phase]["lambda"]
        grids = [read_tsv(run/lang/f"{phase}_lambda_grid.tsv") for lang in outputs]
        assert len(grids[0]) == len(grids[1]) == 100
        for a, b in zip(*grids):
            for column in ("lambda", "absolute_spearman"):
                assert math.isclose(float(a[column]), float(b[column]), rel_tol=1e-10, abs_tol=1e-12)

    differences = [{"source_row": row["source_row"], "gene_id": row["gene_id"],
                    "saved_short_score": row["score_short"],
                    "recomputed_short_score": row["recomputed_score_short"],
                    "recomputed_minus_saved": float(row["recomputed_score_short"])-float(row["score_short"])}
                   for row in outputs["python"]]
    above_precision = [row for row in differences if abs(row["recomputed_minus_saved"]) > 5e-8]
    write_tsv(above_precision, run/"short_recomputed_differences_above_5e-8.tsv", list(differences[0]))
    special = next(row for row in outputs["python"] if row["gene_id"] == "Glyma.14G129500")
    clue = stable_fosco(float(special["raw_short"]), 12, .47)
    report = {"status": "PASS", "records": len(reference),
        "unique_genes": len({row["gene_id"] for row in outputs["python"]}),
        "exact_short_score_text_matches_per_language": len(reference),
        "exact_short_adjusted_text_matches_per_language": len(reference),
        "cross_language_max_absolute_errors": language_errors,
        "recomputed_lambdas": {p: float(summaries["python"][p]["lambda"]) for p in ("short", "middle")},
        "recomputed_short_differences_above_5e_8": len(above_precision),
        "recomputed_short_max_absolute_difference_from_saved": max(abs(r["recomputed_minus_saved"]) for r in differences),
        "special_gene": special,
        "historical_clue": {"adjusted_using_snp12_lambda047": clue,
                            "rounded_to_7_decimals": round(clue, 7),
                            "status": "Numerical consistency only; historical operation unproven"},
        "reference_sha256": sha256(args.reference), "input_sha256": sha256(args.input),
        "source_sha256": provenance["source_sha256"],
        "outputs_sha256": {str(path.relative_to(run)): sha256(path)
                           for lang in outputs for path in sorted((run/lang).glob("*.tsv"))},
        "scope": "Exact preservation of author-designated original short scores, not regeneration of their historical derivation. Middle and diagnostic short recomputed after SNP imputation. No new pathway simulations."}
    destination = ROOT/"reports/paper_short_reference_verification_2026-10-21.json"
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("status", "records", "exact_short_score_text_matches_per_language",
        "cross_language_max_absolute_errors", "recomputed_short_differences_above_5e_8", "special_gene", "historical_clue")}, indent=2))


if __name__ == "__main__":
    main()
