"""Filter full candidate tables using explicit, saved Excel-style rules."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from pathway_selection import select_candidates


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidates", type=Path)
    parser.add_argument("rules", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        parser.error("Output must be new or empty")
    with args.candidates.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rules = json.loads(args.rules.read_text(encoding="utf-8-sig"))
    audit = select_candidates(rows, rules)
    args.output.mkdir(parents=True, exist_ok=True)
    for name, records in (("selection_audit.tsv", audit), ("selected_pathways.tsv", [r for r in audit if r["selected"]])):
        with (args.output / name).open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(audit[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(records)
    counts = {ph: sum(r["selected"] for r in audit if r["phase"] == ph) for ph in ("short", "middle")}
    manifest = dict(status="completed", rules=rules, selected_counts=counts,
                    inputs=[dict(file=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in (args.candidates,args.rules)])
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    print(f"Selected {counts}. Rules: {rules['label']}")


if __name__ == "__main__":
    main()
