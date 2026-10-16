"""Compare R/Python with exactly the same finite null samples, not just a seed."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"python"))
import numpy as np
from ctpathway import read_tsv, monte_carlo_p, maxmean_stat, gsea_es, fosco, fit_fosco


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    out = ROOT/"tests/tmp/cross_language"
    subprocess.run([args.rscript, "--vanilla", "tests/cross_language_reference.R", str(out)],
                   cwd=ROOT, check=True)
    scores = np.array([5, 0, 2, 7, 1, 4, 8, 3, 6, 2], dtype=float)
    hits = [0, 3, 6]  # R positions 1, 4, 7; exactly the same genes.
    null_sets = list(itertools.combinations(range(len(scores)), 3))
    z = (scores-scores.mean())/scores.std(ddof=1)
    ranked = np.sort(scores)[::-1]
    functions = {"sumstat": lambda i: scores[list(i)].sum(),
                 "maxmean": lambda i: maxmean_stat(z[list(i)]),
                 "gsea": lambda i: gsea_es(ranked, i)}
    checks = []
    for row in read_tsv(out/"reference.tsv"):
        statistic = functions[row["method"]]
        observed = statistic(hits)
        p_value = monte_carlo_p(observed, [statistic(i) for i in null_sets])
        np.testing.assert_allclose([observed, p_value], [float(row["observed"]), float(row["p_value"])], rtol=1e-12)
        checks.append({"method": row["method"], "observed": float(observed), "shared_null_p": p_value})
    rows = read_tsv(out/"fosco.tsv")
    raw = [float(r["raw"]) for r in rows]
    snp = [float(r["snp"]) for r in rows]
    np.testing.assert_allclose(fosco(raw, snp, .47), [float(r["adjusted"]) for r in rows], rtol=1e-12)
    fit = fit_fosco(raw, snp)
    np.testing.assert_allclose(fit["adjusted"], [float(r["fitted"]) for r in rows], rtol=1e-12)
    assert np.isclose(fit["lambda"], float(rows[0]["lambda"]))
    report = {"status": "PASS", "null_sets_per_method": len(null_sets), "methods": checks,
              "fosco_fixed_and_fitted": "PASS", "fitted_lambda": fit["lambda"]}
    (ROOT/"reports/cross_language_reference.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print("PASS: R and Python agree on SUMSTAT, GSEA, MaxMean, identical-null p-values, and FOSCO.")


if __name__ == "__main__":
    main()
