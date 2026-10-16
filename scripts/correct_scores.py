"""Correct both phases with FOSCO and save diagnostic plots.

Usage: python scripts/correct_scores.py confirmed_raw_scores.tsv results/python_correction
Input columns: gene_id, raw_short, raw_middle, snp. Raw values must be in (0,1].
This script does not infer probabilities from the collaborator's external scores.
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
from ctpathway import canonical_id, fit_fosco, read_tsv, write_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = read_tsv(args.input)
    required = {"gene_id", "raw_short", "raw_middle", "snp"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Required columns: {sorted(required)}")
    for row in rows:
        row["gene_id"] = canonical_id(row["gene_id"])
    if len({r["gene_id"] for r in rows}) != len(rows) or any(not r["gene_id"] for r in rows):
        raise ValueError("Resolve duplicate or missing gene IDs before correction")
    snp = np.array([float(r["snp"]) for r in rows])
    args.output.mkdir(parents=True, exist_ok=True)
    summaries = []
    for phase in ("short", "middle"):
        raw = np.array([float(r[f"raw_{phase}"]) for r in rows])
        fit = fit_fosco(raw, snp)
        corrected = fit["adjusted"]
        for row, value in zip(rows, corrected):
            row[f"adjusted_{phase}"] = float(value)
            row[f"score_{phase}"] = float(-np.log10(value))
        write_tsv(fit["grid"], args.output/f"{phase}_lambda_grid.tsv")
        summaries.append({"phase": phase, "lambda": fit["lambda"],
                          "raw_spearman": float(spearmanr(snp, raw).statistic),
                          "adjusted_spearman": float(spearmanr(snp, corrected).statistic)})
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
