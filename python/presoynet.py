"""Offline pathway selection, Jaccard crosstalk, and gene-list export.

These functions stop before co-functional network inference. A shared gene
creates a pathway-overlap edge, not evidence of a gene-gene interaction.
"""
from itertools import combinations
import math
from pathlib import Path

from ctpathway import canonical_id, read_gmt, read_tsv, write_tsv

PHASES = ("short", "middle")
METHODS = ("hypergeometric", "sumstat", "gsea")
AUDIT_COLUMNS = ["phase", "go_id", "description", "annotation_size", "tested_size",
                 "ct_count", "ct_ratio", *[m + "_bh" for m in METHODS],
                 "methods_below_alpha", "selected", "selection_basis"]
EDGE_COLUMNS = ["source", "target", "intersection_size", "union_size", "jaccard", "retained", "filter_jaccard"]


def checked_ids(rows, field, label):
    ids = [canonical_id(row[field]) for row in rows]
    if any(not gene for gene in ids) or len(ids) != len(set(ids)):
        raise ValueError(f"Missing or duplicate IDs in {label}")
    return set(ids)


def select_pathways(pathways, tables, universe, foreground, phase, alpha=.05,
                    min_ct=1, min_ratio=0, method_rule="any", requested=None):
    """Select after BH, preserving the original testing family and probabilities.

requested is an explicit set of pathway IDs, or None for a fresh BH selection.
Provided lists bypass the statistical cutoffs and are labelled accordingly.
"""
    if not 0 < alpha <= 1 or min_ct < 0 or int(min_ct) != min_ct or not 0 <= min_ratio <= 1:
        raise ValueError("Invalid selection thresholds")
    if method_rule not in ("any", "all"):
        raise ValueError("method_rule must be any or all")
    by_id = {p["go_id"]: p for p in pathways}
    if len(by_id) != len(pathways):
        raise ValueError("Duplicate GO IDs")
    method_maps = {}
    for method in METHODS:
        rows = tables[method]
        checked_ids(rows, "go_id", method)
        method_maps[method] = {r["go_id"]: r for r in rows}
    tested = set(method_maps[METHODS[0]])
    if any(set(method_maps[m]) != tested for m in METHODS):
        raise ValueError("Methods must contain the same complete tested GO family")
    if not tested <= set(by_id) or (requested is not None and not requested <= set(by_id)):
        raise ValueError("Unknown GO ID in results or supplied selection")
    audit, selected = [], []
    for go_id in sorted(tested | (requested or set())):
        p = by_id[go_id]
        members = set(p["genes"]) & universe
        ct = len(members & foreground)
        ratio = ct / len(members) if members else 0.0
        values = {}
        for method in METHODS:
            row = method_maps[method].get(go_id)
            if row is None:
                values[method] = ""
                continue
            # Reject mixed input versions before using any statistical results.
            for key, expected in (("annotation_size", len(p["genes"])), ("tested_size", len(members)), ("ct_count", ct)):
                if float(row[key]) != expected:
                    raise ValueError(f"Membership mismatch: {phase} {method} {go_id} {key}")
            if not math.isclose(float(row["ct_ratio"]), ratio, rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError(f"CT ratio mismatch: {go_id}")
            value = float(row["hypergeometric_bh" if method == "hypergeometric" else "p_bh"])
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"Invalid BH probability: {method} {go_id}")
            values[method] = value
        passing = [m for m in METHODS if values[m] != "" and values[m] < alpha]
        accepted = (go_id in requested if requested is not None else
                    ct >= min_ct and ratio >= min_ratio and
                    (len(passing) >= (1 if method_rule == "any" else len(METHODS))))
        basis = "provided_list" if requested is not None else "bh_reanalysis"
        audit.append(dict(phase=phase, go_id=go_id, description=p["description"],
                          annotation_size=len(p["genes"]), tested_size=len(members), ct_count=ct,
                          ct_ratio=ratio, **{m+"_bh": values[m] for m in METHODS},
                          methods_below_alpha=";".join(passing), selected=int(accepted), selection_basis=basis))
        if accepted:
            selected.append(p)
    return audit, selected


def jaccard_edges(memberships, keep_fraction):
    """Retain the largest positive undirected overlaps, including cutoff ties.

    Rank only positive pairs; keep ceil(fraction * positive_pairs), then retain
    every pair tied at that cutoff. Self loops and zero-overlap edges are absent.
    """
    if not math.isfinite(keep_fraction) or not 0 <= keep_fraction <= 1:
        raise ValueError("Edge fraction must be between zero and one")
    if any(not genes for genes in memberships.values()):
        raise ValueError("A selected pathway has no genes under this membership policy")
    edges = []
    for a, b in combinations(sorted(memberships), 2):
        common = len(memberships[a] & memberships[b])
        total = len(memberships[a] | memberships[b])
        edges.append(dict(source=a, target=b, intersection_size=common, union_size=total,
                          jaccard=common/total, retained=0))
    positive = sorted((r["jaccard"] for r in edges if r["intersection_size"]), reverse=True)
    count = math.ceil(keep_fraction * len(positive))
    cutoff = positive[count-1] if count else None
    for row in edges:
        row["retained"] = int(cutoff is not None and row["intersection_size"] > 0 and row["jaccard"] >= cutoff)
    return edges, cutoff


def load_phase(analysis, inputs, phase):
    scores = read_tsv(analysis / f"{phase}_analysis_scores.tsv")
    universe = checked_ids(scores, "gene_id", "analysis scores")
    if not universe or any(not math.isfinite(float(r["analysis_score"])) or float(r["analysis_score"]) < 0 for r in scores):
        raise ValueError("Analysis scores must be nonempty and finite")
    foreground = checked_ids(read_tsv(inputs / f"ct_{phase}.tsv"), "gene_id", "CT list")
    tables = {m: read_tsv(analysis / f"{phase}_{m}.tsv") for m in METHODS}
    return universe, foreground, tables


def filter_edges(memberships, keep_fraction, edge_rule, phase):
    """Recalculate coefficients, then apply a stated filtering procedure.

    The dated 20240916 script filters matrices printed to nine decimals using
    fixed cutoffs. Rounding matters at the short-phase boundary (5 / 1122).
    Neither rule uses a target number of retained edges.
    """
    if edge_rule not in ("fraction", "historical-code"):
        raise ValueError("Unknown edge rule")
    edges, cutoff = jaccard_edges(memberships, keep_fraction)
    if edge_rule == "historical-code":
        cutoff = {"short": .004456328, "middle": .14160401}[phase]
    for row in edges:
        row["filter_jaccard"] = round(row["jaccard"], 9) if edge_rule == "historical-code" else row["jaccard"]
        if edge_rule == "historical-code":
            row["retained"] = int(row["filter_jaccard"] >= cutoff)
    return edges, cutoff


def run_presoynet(analysis, inputs, output, *, alpha=.05, min_ct=1, min_ratio=0,
                  method_rule="any", membership="measured", keep_short=.8, keep_middle=.5,
                  selected_pathways=None, selection_label="", edge_rule="fraction"):
    """Calculate both phases before creating outputs, then write auditable tables."""
    analysis, inputs, output = map(Path, (analysis, inputs, output))
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Output must be new or empty")
    if membership not in ("measured", "annotation"):
        raise ValueError("membership must be measured or annotation")
    requested = None
    if selected_pathways is not None:
        if not selection_label.strip():
            raise ValueError("A supplied selection requires --selection-label")
        provided = read_tsv(selected_pathways)
        if any(r["phase"] not in PHASES for r in provided):
            raise ValueError("Selected pathway phase must be short or middle")
        requested = {ph: checked_ids([r for r in provided if r["phase"] == ph], "go_id", "selection") for ph in PHASES}
    pathways = read_gmt(inputs / "pathways.gmt")
    data, summaries = {}, []
    for phase, fraction in zip(PHASES, (keep_short, keep_middle)):
        universe, ct, tables = load_phase(analysis, inputs, phase)
        audit, selected = select_pathways(pathways, tables, universe, ct, phase, alpha, min_ct,
                                         min_ratio, method_rule, None if requested is None else requested[phase])
        members = {p["go_id"]: set(p["genes"]) & universe if membership == "measured" else set(p["genes"]) for p in selected}
        edges, cutoff = filter_edges(members, fraction, edge_rule, phase)
        data[phase] = dict(universe=universe, ct=ct, audit=audit, selected=selected, members=members, edges=edges)
        summaries.append(dict(phase=phase, tested_pathways=len(tables["hypergeometric"]), selected_pathways=len(selected),
                              positive_edges=sum(e["intersection_size"] > 0 for e in edges),
                              retained_edges=sum(e["retained"] for e in edges),
                              keep_fraction=fraction if edge_rule == "fraction" else "", edge_rule=edge_rule,
                              jaccard_cutoff="" if cutoff is None else cutoff))
    output.mkdir(parents=True, exist_ok=True)
    membership_rows = []
    for phase, d in data.items():
        write_tsv(d["audit"], output / f"{phase}_selection_audit.tsv", AUDIT_COLUMNS)
        write_tsv([r for r in d["audit"] if r["selected"]], output / f"{phase}_selected_pathways.tsv", AUDIT_COLUMNS)
        write_tsv(d["edges"], output / f"{phase}_all_pairs.tsv", EDGE_COLUMNS)
        lookup = {(r["source"], r["target"]): r["jaccard"] for r in d["edges"]}
        ids = sorted(d["members"])
        matrix = [dict(go_id=a, **{b: 1.0 if a == b else lookup[tuple(sorted((a, b)))] for b in ids}) for a in ids]
        write_tsv(matrix, output / f"{phase}_jaccard.tsv", ["go_id", *ids])
        membership_rows.extend(dict(phase=phase, go_id=go, gene_id=g, is_ct=int(g in d["ct"]),
                                    is_measured=int(g in d["universe"])) for go in ids for g in sorted(d["members"][go]))
    write_tsv(membership_rows, output / "pathway_genes.tsv", ["phase", "go_id", "gene_id", "is_ct", "is_measured"])
    node_ids = sorted(set(data["short"]["members"]) | set(data["middle"]["members"]))
    annotations = {p["go_id"]: p for p in pathways}
    nodes = []
    for go in node_ids:
        present = [ph for ph in PHASES if go in data[ph]["members"]]
        genes = set().union(*(data[ph]["members"][go] for ph in present))
        nodes.append(dict(go_id=go, description=annotations[go]["description"], phase="both" if len(present) == 2 else present[0],
                          network_gene_count=len(genes), **{ph+"_ct_count": len(data[ph]["members"].get(go, set()) & data[ph]["ct"]) for ph in PHASES}))
    write_tsv(nodes, output / "network_nodes.tsv", ["go_id", "description", "phase", "network_gene_count", "short_ct_count", "middle_ct_count"])
    combined = {}
    for phase, d in data.items():
        for edge in d["edges"]:
            if not edge["retained"]:
                continue
            key = (edge["source"], edge["target"])
            if key not in combined:
                combined[key] = dict(source=key[0], target=key[1], short_jaccard="", middle_jaccard="")
            combined[key][phase+"_jaccard"] = edge["jaccard"]
    write_tsv([combined[k] for k in sorted(combined)], output / "network_edges.tsv", ["source", "target", "short_jaccard", "middle_jaccard"])
    # Gene lists include isolates. These are handoff inputs, not SoyNet outputs.
    export_dir = output / "soynet_inputs"
    export_dir.mkdir()
    shared = set(data["short"]["members"]) & set(data["middle"]["members"])
    groups = {"short": {"short": set(data["short"]["members"])},
              "middle": {"middle": set(data["middle"]["members"])},
              "shared_pathways": {ph: shared for ph in PHASES},
              "combined": {ph: set(data[ph]["members"]) for ph in PHASES}}
    exports = []
    for group, phase_ids in groups.items():
        genes = set().union(*(data[ph]["members"][go] for ph, ids in phase_ids.items() for go in ids))
        ct = genes & set().union(*(data[ph]["ct"] for ph in phase_ids))
        for suffix, values in (("genes", genes), ("ctgenes", ct)):
            (export_dir / f"{group}_{suffix}.txt").write_text("".join(g+"\n" for g in sorted(values)), encoding="utf-8")
        exports.append(dict(group=group, genes=len(genes), ctgenes=len(ct)))
    write_tsv(exports, output / "export_summary.tsv")
    write_tsv(summaries, output / "network_summary.tsv")
    return nodes, list(combined.values()), summaries


def plot_network(nodes, edges, path):
    """Deterministic grouped layout; node placement does not define gene modules."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = {"short": "#377EB8", "both": "#DDAA22", "middle": "#D86999"}
    positions = {}
    for phase, x in (("short", 0), ("both", 1), ("middle", 2)):
        group = sorted(n["go_id"] for n in nodes if n["phase"] == phase)
        for i, go in enumerate(group):
            angle = 2*math.pi*i/len(group)
            positions[go] = ((x+.27*math.cos(angle), .5+.4*math.sin(angle)) if len(group) > 2
                             else (x, 1-(i+1)/(len(group)+1)))
    height = max(6, max([sum(n["phase"] == ph for n in nodes) for ph in colors], default=0)/4)
    fig, ax = plt.subplots(figsize=(12, height))
    for edge in edges:
        a, b = positions[edge["source"]], positions[edge["target"]]
        weight = max(float(edge[k]) for k in ("short_jaccard", "middle_jaccard") if edge[k] != "")
        ax.plot([a[0], b[0]], [a[1], b[1]], color="#777777", alpha=.22, linewidth=.4+2*weight, zorder=1)
    for node in nodes:
        x, y = positions[node["go_id"]]
        ax.scatter(x, y, s=40+2*math.sqrt(node["network_gene_count"]), color=colors[node["phase"]], zorder=2)
        ax.text(x+.035, y, node["go_id"], fontsize=8, va="center")
    for phase, x in (("short", 0), ("both", 1), ("middle", 2)):
        ax.text(x, 1.02, phase.capitalize(), color=colors[phase], ha="center", weight="bold")
    if not nodes:
        ax.text(1, .5, "No pathways pass the selected rule", ha="center")
    ax.set(xlim=(-.5, 2.6), ylim=(-.06, 1.1), title="Pathway crosstalk from gene overlap\nNode area increases with gene count; edges show retained Jaccard overlap")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
