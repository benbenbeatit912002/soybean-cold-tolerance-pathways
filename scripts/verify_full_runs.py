"""Check completed 10,000-draw runs and compare Python/R numerical results.

Run from the repository root after all four full configurations finish:
    python scripts/verify_full_runs.py
Random p-values are compared as Monte Carlo estimates, not for exact equality.
"""
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
import numpy as np
from ctpathway import read_tsv, bh_adjust


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    root = Path(__file__).resolve().parents[1]
    report = {"permutations": 10000, "runs": [], "language_comparisons": [],
              "comparison_note": "Native random generators differ. The six-standard-error screen is diagnostic, not a proof of equivalence.",
              "sha256": {}}
    tables = {}
    for language in ("r", "python"):
        for column in (5, 6):
            name = f"{language}_score{column}_full"
            folder = root/"results"/name
            summary = read_tsv(folder/"input_summary.tsv")  # Written only after both phases finish.
            run = {"name": name, "input_summary": summary, "tests": []}
            for phase, expected_count in (("short", 2638), ("middle", 2637)):
                for method in ("hypergeometric", "sumstat", "gsea", "maxmean"):
                    path = folder/f"{phase}_{method}.tsv"
                    rows = read_tsv(path)
                    tables[language, column, phase, method] = rows
                    check(len(rows) == expected_count, f"{path}: unexpected pathway count")
                    check(len({r['go_id'] for r in rows}) == len(rows), f"{path}: duplicate GO IDs")
                    pcol, qcol = (("hypergeometric_p", "hypergeometric_bh") if method == "hypergeometric"
                                  else ("p_value", "p_bh"))
                    p = np.array([float(r[pcol]) for r in rows])
                    q = np.array([float(r[qcol]) for r in rows])
                    check(np.isfinite(p).all() and np.isfinite(q).all(), f"{path}: nonfinite probabilities")
                    check(np.all((0 <= p) & (p <= q+1e-14) & (q <= 1)), f"{path}: probability bounds")
                    np.testing.assert_allclose(q, bh_adjust(p), rtol=1e-12, atol=1e-14)
                    if method != "hypergeometric":
                        counts = np.array([int(r['null_exceedances']) for r in rows])
                        check(np.all((counts >= 0) & (counts <= 10000)), f"{path}: invalid null counts")
                        check(all(int(r['permutations']) == 10000 for r in rows), f"{path}: incomplete draws")
                        check(all(int(r['seed']) == 20261015+(phase == 'middle') for r in rows), f"{path}: unexpected seed")
                        check(np.isfinite([float(r['statistic']) for r in rows]).all(), f"{path}: invalid statistic")
                        np.testing.assert_allclose(p, (counts+1)/10001, rtol=1e-12, atol=1e-15)
                    run["tests"].append({"phase": phase, "method": method, "pathways": len(rows),
                        "distinct_sizes": len({r['tested_size'] for r in rows}),
                        "min_p": float(p.min()), "min_bh": float(q.min()),
                        "bh_below_0_05": int(np.count_nonzero(q < .05))})
            report["runs"].append(run)
            for path in sorted(folder.iterdir()):
                if path.is_file():
                    report["sha256"][path.relative_to(root).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()

    # Exact counts and deterministic statistics must agree in every real-data row.
    for column in (5, 6):
        r_folder = root/"results"/f"r_score{column}_full"
        py_folder = root/"results"/f"python_score{column}_full"
        check(read_tsv(r_folder/"input_summary.tsv") == read_tsv(py_folder/"input_summary.tsv"),
              f"Input summary differs for score column {column}")
        for phase in ("short", "middle"):
            r_scores = read_tsv(r_folder/f"{phase}_analysis_scores.tsv")
            py_scores = read_tsv(py_folder/f"{phase}_analysis_scores.tsv")
            check([r['gene_id'] for r in r_scores] == [r['gene_id'] for r in py_scores],
                  f"Scored gene universe differs for {phase}, column {column}")
            np.testing.assert_allclose([float(r['analysis_score']) for r in r_scores],
                                       [float(r['analysis_score']) for r in py_scores], rtol=1e-12, atol=1e-14)
            for method in ("hypergeometric", "sumstat", "gsea", "maxmean"):
                r_rows, py_rows = (tables[language, column, phase, method] for language in ("r", "python"))
                for r, p in zip(r_rows, py_rows):
                    check(all(r[key] == p[key] for key in ("go_id", "description", "annotation_size", "tested_size", "ct_count")),
                          f"Mismatched membership: score {column}, {phase}, {method}")
                comparison = {"column": column, "phase": phase, "method": method, "pathways": len(r_rows)}
                if method == "hypergeometric":
                    for key in ("hypergeometric_p", "hypergeometric_bh"):
                        np.testing.assert_allclose([float(r[key]) for r in r_rows],
                                                   [float(p[key]) for p in py_rows], rtol=1e-10, atol=0)
                    comparison["exact_probability_agreement"] = True
                else:
                    r_stat = np.array([float(r['statistic']) for r in r_rows])
                    py_stat = np.array([float(r['statistic']) for r in py_rows])
                    np.testing.assert_allclose(r_stat, py_stat, rtol=1e-10, atol=1e-12)
                    comparison["max_absolute_statistic_difference"] = float(np.max(np.abs(r_stat-py_stat)))
                    r_p = np.array([float(r['p_value']) for r in r_rows])
                    py_p = np.array([float(r['p_value']) for r in py_rows])
                    pooled = (r_p+py_p)/2
                    threshold = 6*np.sqrt(pooled*(1-pooled)*2/10000) + 2/10001
                    comparison["max_absolute_p_difference"] = float(np.max(np.abs(r_p-py_p)))
                    comparison["outside_six_standard_error_screen"] = int(np.count_nonzero(np.abs(r_p-py_p) > threshold))
                report["language_comparisons"].append(comparison)

    for pattern in ("R/*.R", "python/*.py", "scripts/*.R", "scripts/*.py", "config*.R", "config*.json", "requirements.txt"):
        for path in sorted(root.glob(pattern)):
            report["sha256"][path.relative_to(root).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = json.loads((root/"reports/source_manifest.json").read_text(encoding="utf-8"))
    for item in manifest:
        check(hashlib.sha256((root/item['staged']).read_bytes()).hexdigest() == item['staged_sha256'], item['staged'])
    report["verified_staged_hashes"] = len(manifest)
    report["simulation_tables"] = 24
    report["simulation_result_rows"] = sum(len(rows) for key, rows in tables.items() if key[-1] != "hypergeometric")
    report["null_set_evaluations"] = sum(test['distinct_sizes']*10000 for run in report['runs']
                                          for test in run['tests'] if test['method'] != 'hypergeometric')
    report["status"] = "passed"
    path = root/"reports/full_simulation_verification.json"
    path.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(f"PASS: {report['simulation_tables']} simulation tables; {report['simulation_result_rows']:,} pathway results; "
          f"{report['null_set_evaluations']:,} null-set evaluations. See {path}")
    unusual = sum(c.get('outside_six_standard_error_screen', 0) for c in report['language_comparisons'])
    print(f"Cross-language deterministic statistics agree; {unusual} p-value comparisons outside the diagnostic screen.")


if __name__ == "__main__":
    main()
