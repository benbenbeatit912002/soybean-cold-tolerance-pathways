"""Recount primary/measured CTgenes and replay the saved annotation-size test.

This diagnostic reads gene IDs, never saved pathway probabilities. Its historical
hypergeometric calculation uses full annotation sizes and measured CT membership,
as in the archived January script. The main analysis uses measured pathway sizes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from ctpathway import bh_adjust, canonical_id, read_gmt, read_tsv, write_tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene_table", type=Path, help="TSV with gene_id; only this column is used")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        parser.error("Output directory must be new or empty")
    rows = read_tsv(args.gene_table)
    if not rows or any(not r.get("gene_id", "").strip() for r in rows):
        raise ValueError("A nonempty gene_id column is required")
    universe = {canonical_id(r["gene_id"]) for r in rows}
    pathways = read_gmt(ROOT / "data/input/pathways.gmt")
    args.output.mkdir(parents=True, exist_ok=True)
    summaries = []
    missing = []
    sources = [args.gene_table, ROOT / "data/input/pathways.gmt"]
    for phase in ("short", "middle"):
        ct_path = ROOT / f"data/input/ct_{phase}.tsv"
        sources.append(ct_path)
        primary = {canonical_id(r["gene_id"]) for r in read_tsv(ct_path)}
        measured = primary & universe
        missing.extend(dict(phase=phase, gene_id=g) for g in sorted(primary - measured))
        membership = []
        for pathway in pathways:
            genes = set(pathway["genes"])
            membership.append(dict(phase=phase, go_id=pathway["go_id"], description=pathway["description"],
                annotation_size=len(genes), measured_size=len(genes & universe),
                primary_ct_count=len(genes & primary), measured_ct_count=len(genes & measured),
                excluded_ct_ids=";".join(sorted(genes & (primary - measured)))))
        write_tsv(membership, args.output / f"{phase}_membership.tsv")
        candidates = [r for r in membership if 5 <= r["annotation_size"] < 1500 and r["measured_ct_count"] > 0]
        x = np.asarray([r["measured_ct_count"] - 1 for r in candidates])
        sizes = np.asarray([r["annotation_size"] for r in candidates])
        if np.any(sizes > len(universe)):
            raise ValueError("Historical annotation size exceeds the supplied universe")
        stable = hypergeom.sf(x, len(universe), sizes, len(measured))
        # Preserve both expressions so subtraction rounding can be audited.
        legacy = 1 - hypergeom.cdf(x, len(universe), sizes, len(measured))
        calculated = [dict(r, universe_size=len(universe), ct_population=len(measured),
            p_value=float(p), p_bh=float(q), one_minus_cdf=float(old), one_minus_cdf_bh=float(old_q))
            for r, p, q, old, old_q in zip(candidates, stable, bh_adjust(stable), legacy, bh_adjust(legacy))]
        write_tsv(calculated, args.output / f"{phase}_historical_hypergeometric.tsv")
        summaries.append(dict(phase=phase, gene_records=len(rows), unique_genes=len(universe),
            primary_ct=len(primary), measured_ct=len(measured), excluded_ct=len(primary-measured),
            historical_candidates=len(candidates)))
    write_tsv(summaries, args.output / "input_summary.tsv")
    write_tsv(missing, args.output / "missing_ctgenes.tsv", columns=["phase", "gene_id"])
    (args.output / "manifest.json").write_text(json.dumps(dict(
        purpose="Historical membership audit; full annotation sizes are a legacy convention",
        gene_table_fields_used=["gene_id"], expected_results_used_as_inputs=False,
        inputs=[dict(file=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources]
    ), indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summaries))


if __name__ == "__main__":
    main()
