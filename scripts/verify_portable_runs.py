"""Verify two relocated real-data runs using only files inside this repository.

Checks independent hypergeometric/BH calculations, simulation bookkeeping,
short-score reference agreement, cross-language scores and observed statistics.
Native Monte Carlo p-values are reported, not required to be identical.
"""
import argparse
from collections import defaultdict
from decimal import Decimal
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import false_discovery_control, hypergeom

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from ctpathway import read_tsv, read_gmt, canonical_id


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("python_output", type=Path)
    parser.add_argument("r_output", type=Path)
    parser.add_argument("--permutations", type=int, default=10000)
    args = parser.parse_args()
    reference = read_tsv(ROOT / "data/reference/paper_short.tsv")
    pathways = read_gmt(ROOT / "data/input/pathways.gmt")
    report = dict(status="passed", permutations=args.permutations, runs={}, comparisons=[],
                  historical_repairs_confirmed=False, paper_reproduction_verified=False)
    all_tables, all_scores = {}, {}
    for language, folder in (("python", args.python_output), ("r", args.r_output)):
        scores = read_tsv(folder / "scores/recalculated_scores.tsv")
        assert len(scores) == len(reference) == 42575, "Unexpected score record count"
        for actual, expected in zip(scores, reference):
            assert int(actual["source_row"]) == int(expected["source_row"])
            assert actual["gene_id"] == canonical_id(expected["gene_id"])
            assert Decimal(actual["raw_short"]) == Decimal(expected["raw_short"])
            assert Decimal(actual["score_short"]) == Decimal(expected["reference_score_short"])
            assert all(np.isfinite(float(actual[f"score_{phase}"])) for phase in ("short", "middle"))
        all_scores[language] = scores
        run = dict(score_records=len(scores), short_exact_decimal_matches=len(scores), tables=[], random_sets=0)
        for phase in ("short", "middle"):
            by_gene = defaultdict(list)
            for row in scores:
                by_gene[row["gene_id"]].append(float(row[f"score_{phase}"]))
            universe = set(by_gene)
            assert len(universe) == 42567
            means = {gene: float(np.mean(values)) for gene, values in by_gene.items()}
            ct = {r["gene_id"] for r in read_tsv(ROOT / f"data/input/ct_{phase}.tsv")} & universe
            eligible = {p["go_id"]: set(p["genes"]) & universe for p in pathways
                        if 6 <= len(set(p["genes"]) & universe) <= 1500}
            assert len(eligible) == 2689
            methods = ["hypergeometric", "sumstat", "gsea"]
            if (folder / f"analysis/{phase}_maxmean.tsv").exists():
                methods.append("maxmean")
            for method in methods:
                rows = read_tsv(folder / f"analysis/{phase}_{method}.tsv")
                assert len(rows) == len(eligible) and {r["go_id"] for r in rows} == set(eligible)
                for row in rows:
                    genes = eligible[row["go_id"]]
                    assert int(row["tested_size"]) == len(genes)
                    assert int(row["ct_count"]) == len(ct & genes)
                if method == "hypergeometric":
                    p = np.array([float(r["hypergeometric_p"]) for r in rows])
                    q = np.array([float(r["hypergeometric_bh"]) for r in rows])
                    expected = [hypergeom.sf(int(r["ct_count"]) - 1, len(universe), len(ct), int(r["tested_size"])) for r in rows]
                    np.testing.assert_allclose(p, expected, rtol=1e-10, atol=1e-300)
                else:
                    p = np.array([float(r["p_value"]) for r in rows])
                    q = np.array([float(r["p_bh"]) for r in rows])
                    assert all(int(r["permutations"]) == args.permutations for r in rows)
                    exceed = np.array([int(r["null_exceedances"]) for r in rows])
                    assert np.all((exceed >= 0) & (exceed <= args.permutations))
                    np.testing.assert_allclose(p, (exceed + 1) / (args.permutations + 1), rtol=1e-12)
                    if method == "sumstat":
                        expected = [sum(means[g] for g in eligible[r["go_id"]]) for r in rows]
                        np.testing.assert_allclose([float(r["statistic"]) for r in rows], expected, rtol=1e-11, atol=1e-11)
                    run["random_sets"] += len({r["tested_size"] for r in rows}) * args.permutations
                assert np.all(np.isfinite(p)) and np.all((p >= 0) & (p <= 1))
                np.testing.assert_allclose(q, false_discovery_control(p), rtol=1e-10, atol=1e-14)
                assert (folder / f"analysis/{phase}_{method}.png").stat().st_size > 1000
                all_tables[(language, phase, method)] = {r["go_id"]: r for r in rows}
                run["tables"].append(dict(phase=phase, method=method, rows=len(rows)))
        report["runs"][language] = run
    for left, right in zip(all_scores["python"], all_scores["r"]):
        assert left["gene_id"] == right["gene_id"]
        assert left["expression_row"] == right["expression_row"]
        for phase in ("short", "middle"):
            assert Decimal(left[f"score_{phase}"]) == Decimal(right[f"score_{phase}"])
    keys = {(phase, method) for lang, phase, method in all_tables if lang == "python"}
    assert keys == {(phase, method) for lang, phase, method in all_tables if lang == "r"}
    for phase, method in sorted(keys):
        a, b = all_tables[("python", phase, method)], all_tables[("r", phase, method)]
        field = "hypergeometric_p" if method == "hypergeometric" else "statistic"
        np.testing.assert_allclose([float(a[g][field]) for g in a], [float(b[g][field]) for g in a], rtol=1e-10, atol=1e-11)
        item = dict(phase=phase, method=method, max_observed_difference=max(abs(float(a[g][field])-float(b[g][field])) for g in a))
        if method != "hypergeometric":
            item["max_native_monte_carlo_p_difference"] = max(abs(float(a[g]["p_value"])-float(b[g]["p_value"])) for g in a)
        report["comparisons"].append(item)
    report["total_random_sets"] = sum(r["random_sets"] for r in report["runs"].values())
    destination = args.python_output / "portable_verification.json"
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: both recalculated phases, all result tables, BH, simulation counts, and figures. Report: {destination}")


if __name__ == "__main__":
    main()
