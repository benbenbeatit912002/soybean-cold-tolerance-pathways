"""Compare independently generated Python/R crosstalk and historical overlaps."""
import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def compare(left, right):
    a, b = read(left), read(right)
    if len(a) != len(b):
        raise AssertionError(f"Row count differs: {left.name}")
    for x, y in zip(a, b):
        if set(x) != set(y):
            raise AssertionError(f"Columns differ: {left.name}")
        for key in x:
            if x[key] == y[key]:
                continue
            try:
                equal = math.isclose(float(x[key]), float(y[key]), rel_tol=1e-12, abs_tol=1e-12)
            except ValueError:
                equal = False
            if not equal:
                raise AssertionError(f"Values differ: {left.name} {key}: {x[key]} vs {y[key]}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("python_output", type=Path)
    parser.add_argument("r_output", type=Path)
    parser.add_argument("--historical-reference", action="store_true",
                        help="Also compare the historical 13/11 selection and archived Jaccard values")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    files = ["network_nodes.tsv", "network_edges.tsv", "network_summary.tsv", "pathway_genes.tsv", "export_summary.tsv"]
    files += [f"{ph}_{name}.tsv" for ph in ("short", "middle") for name in ("all_pairs", "jaccard")]
    for name in files:
        compare(args.python_output / name, args.r_output / name)
    gene_lists = sorted((args.python_output / "soynet_inputs").glob("*.txt"))
    if len(gene_lists) != 8 or {p.name for p in gene_lists} != {p.name for p in (args.r_output / "soynet_inputs").glob("*.txt")}:
        raise AssertionError("Expected eight matching handoff gene lists")
    for path in gene_lists:
        values = path.read_text(encoding="utf-8").splitlines()
        if values != sorted(set(values)) or values != (args.r_output / "soynet_inputs" / path.name).read_text(encoding="utf-8").splitlines():
            raise AssertionError(f"Unsorted, duplicate or unequal gene export: {path.name}")
    report = dict(status="passed", deterministic_tables_compared=len(files), gene_lists_compared=len(gene_lists),
                  native_simulation_probabilities_compared=False, historical_reference=args.historical_reference)
    if args.historical_reference:
        selected = read(ROOT / "data/reference/crosstalk/selected_pathways.tsv")
        saved = read(ROOT / "data/reference/crosstalk/saved_jaccard_pairs.tsv")
        differences = []
        for phase in ("short", "middle"):
            expected = {r["go_id"] for r in selected if r["phase"] == phase}
            for directory in (args.python_output, args.r_output):
                actual = {r["go_id"] for r in read(directory / f"{phase}_selected_pathways.tsv")}
                if actual != expected:
                    raise AssertionError(f"Historical selection differs: {phase}")
            pairs = {(r["source"], r["target"]): float(r["jaccard"]) for r in read(args.python_output / f"{phase}_all_pairs.tsv")}
            for row in saved:
                if row["phase"] == phase:
                    differences.append(abs(pairs[(row["source"], row["target"])] - float(row["historical_jaccard"])))
        report.update(historical_pairs_compared=len(differences), maximum_jaccard_difference=max(differences),
                      historical_jaccard_matches_saved_precision=max(differences) <= 5.1e-10)
        if not report["historical_jaccard_matches_saved_precision"]:
            raise AssertionError(f"Archived Jaccard values differ: {report}")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
