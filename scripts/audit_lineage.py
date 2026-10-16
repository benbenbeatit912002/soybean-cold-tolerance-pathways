"""Check value relationships between saved artifacts, without claiming biological equivalence."""
import argparse
import json
import math
from collections import defaultdict
from pathlib import Path
from audit_inputs import rows, canonical


def grouped(path, header=False):
    table = rows(path)
    if header:
        next(table)
    result = defaultdict(list)
    for row in table:
        result[canonical(row[0])].append(tuple(float(v) for v in row[1:]))
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("source_root", type=Path)
    args = p.parse_args()
    base = args.source_root / "研究法/大豆基因"
    corrected = grouped(base / "數據/外部基因校正/Yamasaki_Go校正後1213_G.txt", True)
    logged = grouped(base / "數據/外部基因校正/Yamasaki_Go校正後1213_G_log10.txt", True)
    raw = grouped(base / "數據/外部基因校正/gene_p_len.txt", True)
    report = {}
    for label, data in [("legacy_raw", raw), ("legacy_corrected", corrected)]:
        duplicates = {g: v for g, v in data.items() if len(v) > 1}
        report[label] = {"duplicate_ids": len(duplicates),
                         "conflicting_duplicate_ids": sum(len(set(v)) > 1 for v in duplicates.values())}
    for phase, col in [("short", 0), ("middle", 1)]:
        residuals, lambdas = [], []
        for gene in corrected.keys() & logged.keys() & raw.keys():
            if len(corrected[gene]) != 1 or len(logged[gene]) != 1 or len(raw[gene]) != 1:
                continue
            c, l, r = corrected[gene][0], logged[gene][0], raw[gene][0]
            residuals.append(abs(l[col] + math.log10(c[col])))
            if 0 < r[col] < 1 and 0 < c[col] < 1 and c[3] > 1 and r[3] == c[3]:
                lam = math.log(math.log1p(-c[col])/math.log1p(-r[col]))/math.log(c[3])
                if math.isfinite(lam):
                    lambdas.append(lam)
        lambdas.sort()
        ext = grouped(base / f"大豆基因 外部/WITHOUTINF/outer_data_2.0_{phase}_address.txt")
        dups = {g: v for g, v in ext.items() if len(v) > 1}
        report[phase] = {"logged_vs_minus_log10_corrected": {
            "unique_genes_compared": len(residuals),
            "absolute_error_above_1e_minus_6": sum(x > 1e-6 for x in residuals),
            "maximum_absolute_error": max(residuals)},
            "inferred_lambda_from_saved_values": {"n": len(lambdas),
                "median": lambdas[len(lambdas)//2],
                "p05": lambdas[int(len(lambdas)*.05)], "p95": lambdas[int(len(lambdas)*.95)]},
            "external_duplicate_ids": len(dups),
            "external_conflicting_duplicate_ids": sum(len(set(v)) > 1 for v in dups.values()),
            "external_duplicate_examples": dict(list(dups.items())[:3])}
    out = Path(__file__).resolve().parents[1] / "reports/lineage_audit.json"
    out.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
