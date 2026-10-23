"""Verify both real-data correction reruns and compare with the historical file.

Usage: python scripts/verify_snp_correction.py SOURCE_ROOT RUN_DIRECTORY
Run directory contains input.tsv, python/, and r/. Comparisons retain historical
records rather than silently aggregating duplicate genes.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
import numpy as np
from ctpathway import canonical_id, read_gmt, read_tsv, write_tsv


def stable_reference(p, snp, lam):
    return 1.0 if p == 1 else -math.expm1(snp**lam * math.log1p(-p))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("run_directory", type=Path)
    args = parser.parse_args()
    run = args.run_directory
    raw = read_tsv(run/"input.tsv")
    python = read_tsv(run/"python/corrected_scores.tsv")
    r = read_tsv(run/"r/corrected_scores.tsv")
    source_directory = args.source_root/"研究法/大豆基因/數據/外部基因校正"
    source_path = source_directory/"gene_p_len.txt"
    source = read_tsv(source_path)
    historical_path = source_directory/"Yamasaki_Go校正後1213.txt"
    historical = read_tsv(historical_path)
    assert len(raw) == len(python) == len(r) == len(historical) == len(source)
    genes = [canonical_id(row["gene_id"]) for row in raw]
    assert genes == [row["gene_id"] for row in python] == [row["gene_id"] for row in r]
    assert genes == [canonical_id(row["gene_2.0"]) for row in historical]
    assert genes == [canonical_id(row["gene_2.0"]) for row in source]
    assert [int(row["source_row"]) for row in raw] == list(range(2, len(source)+2))
    for column, original in (("raw_short", "p_1"), ("raw_middle", "p_24"),
                             ("snp", "snp"), ("gene_length", "len")):
        assert [float(row[column]) for row in raw] == [float(row[original]) for row in source]
    counts = [float(row["snp"]) for row in raw]
    lengths = [int(row["gene_length"]) for row in raw]
    bins = [(length-1)//2000 for length in lengths]
    positives = [value for value in counts if value > 0]
    expected_counts = counts.copy()
    for i, value in enumerate(counts):
        if value == 0:
            donors = [v for v, group in zip(counts, bins) if v > 0 and group == bins[i]]
            donors = donors or positives
            expected_counts[i] = round(math.fsum(donors)/len(donors))
    for output in (python, r):
        assert [float(row["snp"]) for row in output] == expected_counts
        assert [float(row["snp_original"]) for row in output] == counts
        assert [int(row["gene_length"]) for row in output] == lengths
        assert [int(row["source_row"]) for row in output] == [int(row["source_row"]) for row in raw]
        for phase in ("short", "middle"):
            assert [float(row[f"raw_{phase}"]) for row in output] == [float(row[f"raw_{phase}"]) for row in raw]
    audits = [read_tsv(run/language/"snp_imputation.tsv") for language in ("python", "r")]
    assert len(audits[0]) == len(audits[1]) == counts.count(0)
    for a, b in zip(*audits):
        for key in a:
            if key in ("gene_id", "source"):
                assert a[key] == b[key]
            else:
                assert math.isclose(float(a[key]), float(b[key]), rel_tol=1e-12, abs_tol=1e-12)
    summaries = [{row["phase"]: row for row in read_tsv(run/lang/"correction_summary.tsv")}
                 for lang in ("python", "r")]
    report = {"status": "PASS", "records": len(raw), "unique_genes": len(set(genes)),
              "duplicate_extra_rows_preserved": len(raw)-len(set(genes)),
              "zero_input_rows": counts.count(0), "input_nonzero_mean": math.fsum(positives)/len(positives),
              "imputation": audits[0], "phases": {},
              "raw_source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
              "historical_corrected_sha256": hashlib.sha256(historical_path.read_bytes()).hexdigest(),
              "input_sha256": hashlib.sha256((run/"input.tsv").read_bytes()).hexdigest(),
              "scope": "Real-data SNP imputation and FOSCO rerun. No new pathway permutation p-values."}
    changes = []
    for phase, old_column in (("short", "p_1"), ("middle", "p_24")):
        lam = float(summaries[0][phase]["lambda"])
        assert lam == float(summaries[1][phase]["lambda"])
        grid_python = read_tsv(run/"python"/f"{phase}_lambda_grid.tsv")
        grid_r = read_tsv(run/"r"/f"{phase}_lambda_grid.tsv")
        assert len(grid_python) == len(grid_r) == 100
        for a, b in zip(grid_python, grid_r):
            for key in ("lambda", "absolute_spearman"):
                assert math.isclose(float(a[key]), float(b[key]), rel_tol=1e-10, abs_tol=1e-12)
        language_errors = {}
        for column in (f"adjusted_{phase}", f"score_{phase}"):
            a = np.array([float(row[column]) for row in python])
            b = np.array([float(row[column]) for row in r])
            assert np.isfinite(a).all() and np.isfinite(b).all()
            np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12)
            language_errors[column] = float(np.max(np.abs(a-b)))
        independent = [stable_reference(float(row[f"raw_{phase}"]), snp, lam)
                       for row, snp in zip(raw, expected_counts)]
        adjusted = np.array([float(row[f"adjusted_{phase}"]) for row in python])
        np.testing.assert_allclose(adjusted, independent, rtol=1e-12, atol=0)
        assert np.all((adjusted > 0) & (adjusted <= 1))
        old = np.array([float(row[old_column]) for row in historical])
        errors = np.abs(adjusted-old)
        ordinary = np.array(counts) > 0
        report["phases"][phase] = {
            "lambda": lam, "cross_language_max_absolute_errors": language_errors,
            "independent_formula_max_absolute_error": float(np.max(np.abs(adjusted-independent))),
            "historical_corrected_rows_with_absolute_difference_above_1e_12": int(sum(errors > 1e-12)),
            "historical_positive_snp_rows_max_absolute_difference": float(max(errors[ordinary])),
            "all_new_scores_finite": True, "historical_zero_adjusted_rows": int(sum(old == 0))}
        for i in range(len(raw)):
            if counts[i] == 0 or old[i] == 0 or errors[i] > 1e-12:
                previous_snp = float(historical[i]["snp"])
                model_old_snp = stable_reference(float(raw[i][f"raw_{phase}"]), previous_snp, lam)
                changes.append({"source_row": raw[i]["source_row"], "gene_id": genes[i], "phase": phase,
                    "input_snp": counts[i], "new_snp": expected_counts[i], "historical_snp": previous_snp,
                    "lambda": lam, "new_adjusted": adjusted[i], "historical_adjusted": old[i],
                    "new_score": float(python[i][f"score_{phase}"]),
                    "historical_score_from_stored_probability": -math.log10(old[i]) if old[i] > 0 else "Inf",
                    "model_score_using_historical_snp": -math.log10(model_old_snp),
                    "score_change_due_only_to_snp_at_same_lambda":
                        float(python[i][f"score_{phase}"]) + math.log10(model_old_snp)})
    write_tsv(changes, run/"historical_gene_comparison.tsv")
    # Membership-only hypergeometric input is unchanged by replacing an SNP count.
    report["gene_universe_unchanged"] = True
    imputed_genes = {genes[i] for i, value in enumerate(counts) if value == 0}
    terms = read_gmt(Path(__file__).resolve().parents[1]/"data/input/pathways.gmt")
    report["go_terms_containing_imputed_genes"] = [
        p["go_id"] for p in terms if imputed_genes.intersection(p["genes"])]
    report["changed_gene_details"] = changes
    report["outputs_sha256"] = {
        str(path.relative_to(run)): hashlib.sha256(path.read_bytes()).hexdigest()
        for language in ("python", "r") for path in sorted((run/language).glob("*.tsv"))}
    destination = Path(__file__).resolve().parents[1]/"reports/snp_imputation_verification_2026-10-20.json"
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({key: report[key] for key in ("status", "records", "unique_genes", "imputation", "phases",
                                                 "go_terms_containing_imputed_genes", "changed_gene_details")},
                     indent=2))


if __name__ == "__main__":
    main()
