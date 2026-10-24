"""Select pathways, calculate crosstalk, and export gene lists entirely offline."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from presoynet import PHASES, METHODS, run_presoynet, plot_network


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis", type=Path, help="Existing analysis directory with both phases and all three main methods")
    parser.add_argument("output", type=Path, help="New or empty output directory")
    parser.add_argument("--input-dir", type=Path, default=ROOT / "data/input")
    parser.add_argument("--alpha", type=float, default=.05)
    parser.add_argument("--min-ct", type=int, default=1)
    parser.add_argument("--min-ratio", type=float, default=0)
    parser.add_argument("--method-rule", choices=("any", "all"), default="any")
    parser.add_argument("--membership", choices=("measured", "annotation"), default="measured")
    parser.add_argument("--keep-short", type=float, default=.8)
    parser.add_argument("--keep-middle", type=float, default=.5)
    parser.add_argument("--selected-pathways", type=Path, help="Optional phase/go_id TSV; bypasses statistical selection")
    parser.add_argument("--selection-label", default="", help="Required provenance label for a supplied pathway list")
    args = parser.parse_args()
    sources = [args.input_dir / "pathways.gmt"]
    for phase in PHASES:
        sources += [args.input_dir / f"ct_{phase}.tsv", args.analysis / f"{phase}_analysis_scores.tsv"]
        sources += [args.analysis / f"{phase}_{method}.tsv" for method in METHODS]
    if args.selected_pathways:
        sources.append(args.selected_pathways)
    fingerprints = [dict(file=str(p.resolve()), md5=hashlib.md5(p.read_bytes()).hexdigest(),
                         sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources]
    nodes, edges, summary = run_presoynet(args.analysis, args.input_dir, args.output,
        **{k: getattr(args, k) for k in ("alpha", "min_ct", "min_ratio", "method_rule", "membership",
                                        "keep_short", "keep_middle", "selected_pathways", "selection_label")})
    plot_network(nodes, edges, args.output / "crosstalk.png")
    manifest = dict(status="completed", language="Python", selection_mode="provided_list" if args.selected_pathways else "bh_reanalysis",
                    settings={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                    inputs=fingerprints, network_summary=summary, paper_reproduction_verified=False,
                    soynet_contacted=False, edge_rule="ceil(fraction * positive undirected pairs), include all cutoff ties")
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    print(f"Completed offline crosstalk: {len(nodes)} nodes, {len(edges)} edges. Gene lists: {args.output / 'soynet_inputs'}")


if __name__ == "__main__":
    main()
