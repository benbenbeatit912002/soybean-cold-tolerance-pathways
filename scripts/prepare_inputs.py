"""Stage explicitly selected inputs with hashes; do not infer score semantics."""
import argparse
import csv
import hashlib
import json
import shutil
import re
from collections import defaultdict, Counter
from pathlib import Path
from audit_inputs import rows, canonical


def resolve_gene(original, mapping):
    original = original.strip()
    gene = canonical(original)
    if not original:
        return "", "excluded_blank_id"
    if re.fullmatch(r"Glyma\.(?:\d{2}G|u)\d{6}", gene):
        return gene, "v2_id"
    candidates = mapping.get(original, set())
    if len(candidates) == 1:
        return next(iter(candidates)), "mapped_v1_to_v2"
    return "", "excluded_ambiguous_mapping" if candidates else "excluded_unmapped_id"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_root", type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    base = args.source_root / "研究法/大豆基因"
    dest = root / "data/input"
    dest.mkdir(parents=True, exist_ok=True)
    manifest = []
    correspondence = base / "大豆基因 外部/Glyma_genes_correspondence_final.txt"
    mapping = defaultdict(set)
    table = rows(correspondence)
    next(table)
    for row in table:
        target_id = canonical(row[2])
        if re.fullmatch(r"Glyma\.(?:\d{2}G|u)\d{6}", target_id):
            for old_id in row[:2]:
                if old_id and old_id != "no correspondence":
                    mapping[old_id].add(target_id)
    import_summary = {}

    def record(source, target, transform):
        manifest.append({"source": source.relative_to(args.source_root).as_posix(),
                         "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                         "staged": target.relative_to(root).as_posix(),
                         "staged_sha256": hashlib.sha256(target.read_bytes()).hexdigest(),
                         "transformation": transform})

    target = dest / "gene_correspondence.tsv"
    shutil.copyfile(correspondence, target)
    record(correspondence, target, "Unmodified genome-version correspondence table")

    for phase in ("short", "middle"):
        source = base / f"CTgenes 含分數/{phase.upper()}.txt"
        target = dest / f"ct_{phase}.tsv"
        with target.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id", "prior_score"])
            writer.writerows((canonical(r[0]), r[1]) for r in rows(source))
        record(source, target, "Header added; Glyma.nn g/G normalized to G; original prior scores retained")
        source = base / f"大豆基因 外部/WITHOUTINF/outer_data_2.0_{phase}_address.txt"
        target = dest / f"external_{phase}.tsv"
        imported, decisions = [], []
        for source_row, r in enumerate(rows(source), start=1):
            original = r[0].strip()
            gene, status = resolve_gene(original, mapping)
            decisions.append([source_row, original, gene, status])
            if gene:
                imported.append([gene, *r[1:]])
        with target.open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene_id"] + [f"source_col_{c}" for c in range(2, 7)])
            writer.writerows(imported)
        record(source, target, "Neutral headers; v1 IDs mapped only with a unique v2 target; blank/ambiguous/unmapped excluded and audited; duplicate rows retained")
        with (root / f"reports/{phase}_import_decisions.tsv").open("w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(["source_row", "original_gene_id", "gene_id", "status"])
            writer.writerows(decisions)
        import_summary[phase] = {"raw_rows": len(decisions), "retained_rows": len(imported),
                                 "unique_gene_ids": len({r[0] for r in imported}),
                                 "statuses": dict(Counter(r[3] for r in decisions))}

    source = base / "數據/GO_terms_4899.txt"
    target = dest / "pathways.gmt"
    with target.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        table = rows(source)
        next(table)
        for r in table:
            writer.writerow([r[0], r[1]] + [canonical(g) for g in r[3:] if g.strip()])
    record(source, target, "Wide padded GO table converted to GMT; declared counts verified in audit; gene IDs normalized")

    legacy = root / "legacy"
    legacy.mkdir(exist_ok=True)
    for source in sorted((base / "程式").iterdir()):
        if source.is_file():
            target = legacy / source.name
            shutil.copyfile(source, target)
            record(source, target, "Unmodified archival copy; not executed by the cleaned workflow")
    (root / "reports/source_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    (root / "reports/import_summary.json").write_text(
        json.dumps(import_summary, indent=2)+"\n", encoding="utf-8")
    print(f"Staged {len(manifest)} inputs and archived scripts; source files unchanged.")


if __name__ == "__main__":
    main()
