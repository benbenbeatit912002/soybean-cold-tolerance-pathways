"""Replay identifiable historical SUMSTAT and hypergeometric settings.

This is a forensic comparison, not the recommended analysis. SUMSTAT retains
duplicate records, selects full-annotation sizes 5..1499 with a CT overlap,
and uses strict exceedances / B. Hypergeometric variants expose the historical
matrix-length/size-column mistakes instead of silently adopting them.
No paper results or archived gene scores are used as calculation inputs.
"""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
import numpy as np
from scipy.stats import hypergeom
from ctpathway import read_tsv, read_gmt, write_tsv, bh_adjust


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scores", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261021)
    args = parser.parse_args()
    if args.permutations < 1:
        raise ValueError("permutations must be positive")
    rows = read_tsv(args.scores)
    genes = [r["gene_id"] for r in rows]
    universe = set(genes)
    pathways = read_gmt("data/input/pathways.gmt")
    args.output.mkdir(parents=True, exist_ok=True)
    summaries = []
    for phase in ("short", "middle"):
        ct = read_tsv(f"data/input/ct_{phase}.tsv")
        foreground = {r["gene_id"] for r in ct}
        candidates = [p for p in pathways if 5 <= len(p["genes"]) < 1500 and set(p["genes"]) & foreground]
        # This reproduces length(unique(answer[,3:ncol(answer)])): sizes
        # were counted as if they were gene IDs. It is intentionally labeled.
        annotation_union = set().union(*(set(p["genes"]) for p in candidates))
        mixed_union = annotation_union | {str(len(p["genes"])) for p in candidates}
        values = np.asarray([float(r[f"score_{phase}"]) for r in rows])
        if not np.isfinite(values).all():
            raise ValueError("Historical diagnostics require finite recalculated scores")
        memberships = [set(p["genes"]) for p in candidates]
        indices = [np.asarray([i for i, g in enumerate(genes) if g in members]) for members in memberships]
        sizes = np.asarray([len(index) for index in indices])
        observed = np.asarray([values[index].sum() for index in indices])
        counts = np.empty(len(candidates), dtype=int)
        seed = args.seed + (phase == "middle")
        rng = np.random.default_rng(seed)
        for number, size in enumerate(sorted(set(sizes)), start=1):
            null = np.asarray([values[rng.choice(len(values), size=int(size), replace=False)].sum()
                               for _ in range(args.permutations)])
            for index in np.flatnonzero(sizes == size):
                counts[index] = np.count_nonzero(null > observed[index])
            if number % 25 == 0:
                print(f"{phase}: historical SUMSTAT {number}/{len(set(sizes))} sizes", flush=True)
        pvalues = counts / args.permutations
        base = [dict(go_id=p["go_id"], description=p["description"], annotation_size=len(p["genes"]),
                     tested_size=int(size), ct_count=len(members & foreground))
                for p, size, members in zip(candidates, sizes, memberships)]
        write_tsv([dict(r, statistic=float(statistic), null_exceedances=int(count), p_value=float(p), p_bh=float(q),
                        permutations=args.permutations, seed=seed)
                   for r, statistic, count, p, q in zip(base, observed, counts, pvalues, bh_adjust(pvalues))],
                  args.output/f"{phase}_sumstat.tsv")
        profiles = {"legacy_container_bug": (len(mixed_union), 2*len(foreground)),
                    "genome_with_primary_ct": (len(universe), len(foreground)),
                    "genome_with_measured_ct": (len(universe), len(universe & foreground))}
        hyper_rows = []
        for name, (N, M) in profiles.items():
            probability = np.asarray([hypergeom.sf(r["ct_count"]-1, N, r["annotation_size"], M) for r in base])
            # Keep a stable tail and the cancellation-prone legacy expression.
            legacy_p = np.asarray([1-hypergeom.cdf(r["ct_count"]-1, N, r["annotation_size"], M) for r in base])
            hyper_rows.extend(dict(r, profile=name, universe_size=N, ct_population=M, p_value=float(p),
                                   p_bh=float(q), one_minus_cdf=float(old))
                              for r, p, q, old in zip(base, probability, bh_adjust(probability), legacy_p))
        write_tsv(hyper_rows, args.output/f"{phase}_hypergeometric_variants.tsv")
        summaries.append(dict(phase=phase, records=len(rows), unique_genes=len(universe),
                              candidates=len(candidates), distinct_record_sizes=len(set(sizes)),
                              annotation_union=len(annotation_union), mixed_union=len(mixed_union),
                              ct_input=len(foreground), ct_measured=len(universe & foreground)))
    write_tsv(summaries, args.output/"input_summary.tsv")
    (args.output/"settings.json").write_text(json.dumps(dict(permutations=args.permutations, seed=args.seed,
        score_input=str(args.scores), method="Historical SUMSTAT replay and explicitly invalid hypergeometric diagnostic profiles",
        paper_result_inputs=False, duplicate_policy="retain records", strict_exceedance=True,
        no_exceedance_pseudocount=True), indent=2)+"\n")
    print(f"Outputs: {args.output}")


if __name__ == "__main__":
    main()
