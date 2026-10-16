"""Plot method results: python scripts/plot_results.py results/python_score5_full.

These plots show the strongest 15 terms containing CTgenes, whether significant or not.
"""
import argparse
import os
from pathlib import Path
import sys
import tempfile
import textwrap

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir())/"soybean_matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ctpathway import read_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path, default=Path("results/python_external_membership"))
    args = parser.parse_args()
    colors = {"hypergeometric": "#337A8C", "sumstat": "#527A4D", "gsea": "#76558C", "maxmean": "#B47043"}
    count = 0
    for phase in ("short", "middle"):
        for method, color in colors.items():
            source = args.output/f"{phase}_{method}.tsv"
            if not source.exists():
                continue
            rows = [row for row in read_tsv(source) if int(row["ct_count"]) > 0]
            p_key, q_key = ("hypergeometric_p", "hypergeometric_bh") if method == "hypergeometric" else ("p_value", "p_bh")
            rows = sorted(rows, key=lambda r: (float(r[q_key]), float(r[p_key]), r["go_id"]))[:15]
            if not rows:
                continue
            rows.reverse()
            values = [-np.log10(max(float(r[q_key]), np.finfo(float).tiny)) for r in rows]
            labels = [r["go_id"] + "  " + textwrap.shorten(r["description"], width=49, placeholder="...") for r in rows]
            figure, axis = plt.subplots(figsize=(11, 7))
            figure.subplots_adjust(left=.47, right=.97, top=.9, bottom=.13)
            axis.barh(range(len(rows)), values, color=color)
            axis.set_yticks(range(len(rows)), labels, fontsize=8)
            axis.set_xlabel("-log10(BH-adjusted p-value)")
            title = f"{'Short-term' if phase == 'short' else 'Mid-term'} | {method.upper()}"
            if method == "maxmean":
                title += " (experimental)"
            axis.set_title(title, fontsize=14, weight="bold", pad=15)
            axis.spines[["top", "right"]].set_visible(False)
            figure.text(.5, .025, "Top 15 terms containing CTgenes; selection does not imply significance.",
                        ha="center", fontsize=8)
            figure.savefig(args.output/f"{phase}_{method}.png", dpi=150)
            plt.close(figure)
            count += 1
    if count == 0:
        raise ValueError("No results to plot: run the analysis first and choose its output folder")
    print(f"Saved {count} plots to {args.output}")


if __name__ == "__main__":
    main()
