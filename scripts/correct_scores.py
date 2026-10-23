"""Correct both phases with FOSCO and save diagnostic plots.

Usage: python scripts/correct_scores.py confirmed_raw_scores.tsv results/python_correction
Input columns: gene_id, raw_short, raw_middle, snp. Raw values must be in (0,1].
Include gene_length (bp) when SNP counts contain zero. Zeros use the same-group
nonzero mean, or the global nonzero mean if the group has no donors.
This script does not infer probabilities from the collaborator's external scores.
For historical arithmetic, use reproduce_historical_short.py. This command
always calculates its output scores and never copies a saved score reference.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir())/"soybean_matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import spearmanr
from ctpathway import canonical_id, fit_fosco, impute_snp_counts, read_tsv, write_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--duplicate-policy", choices=("error", "preserve"), default="error",
                        help="preserve retains historical records for a correction audit")
    args = parser.parse_args()
    rows = read_tsv(args.input)
    required = {"gene_id", "raw_short", "raw_middle", "snp"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Required columns: {sorted(required)}")
    for row in rows:
        row["gene_id"] = canonical_id(row["gene_id"])
    if any(not r["gene_id"] for r in rows):
        raise ValueError("Missing gene IDs")
    unique_genes = len({r["gene_id"] for r in rows})
    if unique_genes != len(rows) and args.duplicate_policy == "error":
        raise ValueError("Resolve duplicate gene IDs or explicitly use --duplicate-policy preserve")
    original_snp = np.array([float(r["snp"]) for r in rows])
    lengths = [float(r["gene_length"]) for r in rows] if "gene_length" in rows[0] else None
    imputation = impute_snp_counts(original_snp, lengths)
    snp = imputation["values"]
    for index, row in enumerate(rows):
        row["snp_original"] = float(original_snp[index])
        row["snp"] = float(snp[index])
        row["snp_imputation_source"] = "unchanged"
    for entry in imputation["audit"]:
        row = rows[entry["row_number"] - 1]
        entry["gene_id"] = row["gene_id"]
        row["snp_imputation_source"] = entry["source"]
    args.output.mkdir(parents=True, exist_ok=True)
    write_tsv(imputation["audit"], args.output/"snp_imputation.tsv",
              ["row_number", "gene_id", "bin_start", "bin_end_exclusive", "original_snp",
               "imputed_snp", "donor_count", "donor_mean", "source"])
    settings = {"input": str(args.input), "rows": len(rows), "unique_genes": unique_genes,
                "duplicate_policy": args.duplicate_policy, "bin_width": 2000, "bin_origin": 1,
                "zero_policy": "group_nonzero_mean_then_global_nonzero_mean",
                "donors": "positive_input_rows_before_any_imputation",
                "rounding": "nearest_integer_ties_to_even", "imputed_rows": len(imputation["audit"])}
    settings.update(short_score_source="recomputed_fosco")
    write_tsv([{"setting": key, "value": value} for key, value in settings.items()],
              args.output/"correction_settings.tsv")
    summaries = []
    for phase in ("short", "middle"):
        raw = np.array([float(r[f"raw_{phase}"]) for r in rows])
        fit = fit_fosco(raw, snp)
        corrected = fit["adjusted"]
        for row, value in zip(rows, corrected):
            row[f"adjusted_{phase}"] = float(value)
            row[f"score_{phase}"] = float(-np.log10(value))
        score_source = "recomputed_fosco"
        if phase == "short":
            for row in rows:
                row["short_score_source"] = score_source
        write_tsv(fit["grid"], args.output/f"{phase}_lambda_grid.tsv")
        summaries.append({"phase": phase, "lambda": fit["lambda"],
                          "score_source": score_source,
                          "raw_spearman": float(spearmanr(snp, raw).statistic),
                          "adjusted_spearman": float(spearmanr(snp,
                              [float(row[f"adjusted_{phase}"]) for row in rows]).statistic),
                          "recomputed_adjusted_spearman": float(spearmanr(snp, corrected).statistic)})
        figure, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
        for column, values, label, color in [(0, raw, "raw", "#337A8C"),
                                               (1, corrected, "adjusted", "#CB7551")]:
            negative_log = -np.log10(values)
            axes[0, column].hist(negative_log, bins=min(50, max(5, len(rows)//4)), color=color)
            axes[0, column].set(xlabel=f"-log10({label} score)", ylabel="Frequency",
                                title=f"{phase.capitalize()}: {label} score")
            axes[1, column].scatter(snp, negative_log, s=18 if len(rows) <= 1000 else 3,
                                     alpha=.6, color=color)
            axes[1, column].set(xlabel="SNP count", ylabel=f"-log10({label} score)")
        figure.suptitle(f"FOSCO diagnostics | fitted lambda = {fit['lambda']:.2f}")
        figure.savefig(args.output/f"{phase}_correction_diagnostics.png", dpi=150)
        plt.close(figure)
    write_tsv(rows, args.output/"corrected_scores.tsv")
    write_tsv(summaries, args.output/"correction_summary.tsv")
    (args.output/"environment.json").write_text(json.dumps({"python": sys.version,
        "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__}, indent=2)+"\n")
    print(f"Both phases corrected: {args.output.resolve()}")


if __name__ == "__main__":
    main()
