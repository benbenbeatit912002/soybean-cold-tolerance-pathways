"""Inventory source inputs without modifying them; Python standard library only."""
import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        yield from csv.reader(handle, delimiter="\t")


def canonical(value):
    return re.sub(r"^Glyma\.(\d{2})g", r"Glyma.\1G", value.strip(), flags=re.I)


def profile(path, header=False):
    content = list(rows(path))
    names = content.pop(0) if header else None
    ids = [canonical(r[0]) for r in content]
    stats = {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
             "rows": len(content), "columns": dict(Counter(map(len, content))),
             "unique_ids": len(set(ids)), "duplicate_rows_by_id": len(ids)-len(set(ids)),
             "blank_id_rows": sum(not g for g in ids),
             "unique_nonblank_ids": len({g for g in ids if g}),
             "non_v2_nonblank_id_rows": sum(bool(g) and not bool(re.fullmatch(r"Glyma\.(?:\d{2}G|u)\d{6}", g)) for g in ids),
             "header": names, "numeric_columns": {}}
    for col in range(1, min(map(len, content))):
        values = []
        invalid = 0
        for row in content:
            try:
                values.append(float(row[col]))
            except ValueError:
                invalid += 1
        finite = [v for v in values if math.isfinite(v)]
        stats["numeric_columns"][str(col+1)] = dict(
            minimum=min(finite) if finite else None, maximum=max(finite) if finite else None,
            nonfinite=len(values)-len(finite), nonnumeric=invalid,
            zero=sum(v == 0 for v in finite), negative=sum(v < 0 for v in finite),
            above_one=sum(v > 1 for v in finite))
    return stats, content


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_root", type=Path, help="Original project root")
    parser.add_argument("--output", type=Path, default=Path("reports/input_audit.json"))
    args = parser.parse_args()
    base = args.source_root / "研究法/大豆基因"
    paths = {
        "ct_short": (base / "CTgenes 含分數/SHORT.txt", False),
        "ct_middle": (base / "CTgenes 含分數/MIDDLE.txt", False),
        "external_short": (base / "大豆基因 外部/WITHOUTINF/outer_data_2.0_short_address.txt", False),
        "external_middle": (base / "大豆基因 外部/WITHOUTINF/outer_data_2.0_middle_address.txt", False),
        "legacy_corrected_log10": (base / "數據/外部基因校正/Yamasaki_Go校正後1213_G_log10.txt", True),
        "legacy_raw": (base / "數據/外部基因校正/gene_p_len.txt", True),
        "legacy_corrected": (base / "數據/外部基因校正/Yamasaki_Go校正後1213_G.txt", True),
        "legacy_ct_short": (base / "數據/CT_genes/short_CT_G_score.txt", False),
        "legacy_ct_middle": (base / "數據/CT_genes/middle_CT_G_score.txt", False),
    }
    report, tables = {}, {}
    for key, (path, header) in paths.items():
        report[key], tables[key] = profile(path, header)
        report[key]["path"] = str(path.relative_to(args.source_root)).replace("\\", "/")
    gene_sets = {k: {canonical(r[0]) for r in v} for k, v in tables.items()}
    report["overlaps"] = {}
    for phase in ("short", "middle"):
        ct, ext = gene_sets["ct_"+phase], gene_sets["external_"+phase]
        report["overlaps"][phase] = {
            "ct_in_external": len(ct & ext), "ct_missing_external": sorted(ct-ext),
            "ct_missing_legacy_corrected": sorted(ct-gene_sets["legacy_corrected_log10"]),
            "legacy_ct_id_difference": sorted(ct ^ gene_sets["legacy_ct_"+phase]),
            "external_vs_legacy_common": len(ext & gene_sets["legacy_corrected_log10"]),
        }
    report["overlaps"]["ct_union"] = len(gene_sets["ct_short"] | gene_sets["ct_middle"])
    report["overlaps"]["external_common"] = len(gene_sets["external_short"] & gene_sets["external_middle"])
    go_path = base / "數據/GO_terms_4899.txt"
    go_rows = list(rows(go_path))
    first = go_rows[0][:5]
    if not go_rows[0][0].startswith("GO:"):
        go_rows.pop(0)
    pathways = [(r[0], int(r[2]), {canonical(g) for g in r[3:] if g.strip()}) for r in go_rows]
    report["go"] = {
        "sha256": hashlib.sha256(go_path.read_bytes()).hexdigest(), "first_fields": first,
        "rows": len(pathways), "unique_go_ids": len({p[0] for p in pathways}),
        "unique_genes": len(set().union(*(p[2] for p in pathways))),
        "declared_size_mismatches": sum(n != len(gs) for _, n, gs in pathways),
        "historical_candidates": {phase: sum(4 < n < 1500 and bool(gs & gene_sets["ct_"+phase]) for _, n, gs in pathways) for phase in ("short", "middle")},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
