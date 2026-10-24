"""Recalculate the differences between the saved network procedures.

This audit consumes packaged annotations and historical evidence. No output
count is used as a selection target, and the current 11-pathway list is intact.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"python"))
from ctpathway import read_gmt, read_tsv, write_tsv
from presoynet import filter_edges


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        parser.error("Output must be new or empty")
    reference=ROOT/"data/reference/crosstalk"
    pathways={p["go_id"]:set(p["genes"]) for p in read_gmt(ROOT/"data/input/pathways.gmt")}
    selected=read_tsv(reference/"selected_pathways.tsv")
    differences=[]; summaries=[]; retained={"fraction":set(),"historical-code":set()}
    for phase,fraction in (("short",.8),("middle",.5)):
        members={r["go_id"]:pathways[r["go_id"]] for r in selected if r["phase"]==phase}
        fresh, fresh_cutoff=filter_edges(members,fraction,"fraction",phase)
        old, old_cutoff=filter_edges(members,fraction,"historical-code",phase)
        for a,b in zip(fresh,old):
            key=(a["source"],a["target"])
            for mode,row in (("fraction",a),("historical-code",b)):
                if row["retained"]: retained[mode].add(key)
            if a["retained"] != b["retained"]:
                differences.append(dict(phase=phase,source=key[0],target=key[1],
                    jaccard=a["jaccard"],historical_filter_jaccard=b["filter_jaccard"],
                    fraction_retained=a["retained"],historical_retained=b["retained"]))
        summaries.append(dict(phase=phase,positive_edges=sum(r["jaccard"]>0 for r in fresh),
            fraction_edges=sum(r["retained"] for r in fresh),historical_edges=sum(r["retained"] for r in old),
            fraction_cutoff=fresh_cutoff,historical_cutoff=old_cutoff))
    old_ids={r["go_id"] for r in read_tsv(reference/"historical_middle_14.tsv")}
    current_ids={r["go_id"] for r in selected if r["phase"]=="middle"}
    ct={r["gene_id"] for r in read_tsv(ROOT/"data/input/ct_middle.tsv")}
    old_genes=set().union(*(pathways[g] for g in old_ids)) & ct
    current_genes=set().union(*(pathways[g] for g in current_ids)) & ct
    saved={r["gene_id"] for r in read_tsv(reference/"saved_middle_122_ctgenes.tsv")}
    if old_genes != saved:
        raise AssertionError("Recomputed historical middle membership differs from the independently saved CT list")
    gene_rows=[dict(gene_id=g,in_current=int(g in current_genes),in_historical=int(g in old_genes),
        historical_supporting_go_ids=";".join(sorted(go for go in old_ids if g in pathways[go])),
        current_supporting_go_ids=";".join(sorted(go for go in current_ids if g in pathways[go])))
        for g in sorted(old_genes | current_genes)]
    report=dict(status="passed",fraction_merged_edges=len(retained["fraction"]),
        historical_merged_edges=len(retained["historical-code"]),changed_edges=len(differences),
        current_middle_pathways=len(current_ids),historical_middle_pathways=len(old_ids),
        current_middle_ctgenes=len(current_genes),historical_middle_ctgenes=len(old_genes),
        historical_gene_set_matches_saved=True,only_in_historical=sorted(old_genes-current_genes),
        only_in_current=sorted(current_genes-old_genes),removed_middle_pathways=sorted(old_ids-current_ids))
    args.output.mkdir(parents=True,exist_ok=True)
    write_tsv(differences,args.output/"edge_differences.tsv")
    write_tsv(summaries,args.output/"edge_summary.tsv")
    write_tsv(gene_rows,args.output/"middle_gene_lineage.tsv")
    (args.output/"trace_summary.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__": main()
