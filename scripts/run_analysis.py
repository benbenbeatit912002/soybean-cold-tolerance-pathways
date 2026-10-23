"""Python entry point: load inputs -> audit membership -> test pathways -> save.

From the repository root:
    python scripts/run_analysis.py
    python scripts/run_analysis.py config.score5-full.json
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import platform
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
import numpy as np
import scipy
from ctpathway import (read_tsv, write_tsv, read_gmt, prepare_membership,
                      hypergeometric, bh_adjust, transform_score,
                      resolve_duplicate_scores, permutation_test)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", nargs="?", default="config.python.json")
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    inputs, output = Path(config["input_dir"]), Path(config["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    (output/"settings.json").write_text(json.dumps(config, indent=2)+"\n", encoding="utf-8")
    environment = {"python": sys.version, "platform": platform.platform(),
                   "numpy": np.__version__, "scipy": scipy.__version__,
                   "random_generator": "numpy.random.default_rng (PCG64)"}
    (output/"environment.json").write_text(json.dumps(environment, indent=2)+"\n", encoding="utf-8")
    pathways = read_gmt(inputs/"pathways.gmt")
    summaries = []
    for phase in config["phases"]:
        if phase not in ("short", "middle"):
            raise ValueError("Only short and middle phases are supported")
        # 1. Keep CT membership distinct from its prior numerical score.
        ct = read_tsv(inputs/f"ct_{phase}.tsv")
        score_input = config.get("score_input")
        external = read_tsv(Path(score_input) if score_input else inputs/f"external_{phase}.tsv")
        if not ct or not external or any(not r["gene_id"] for r in ct+external):
            raise ValueError("Missing gene identifiers")
        foreground = [row["gene_id"] for row in ct]
        if len(set(foreground)) != len(foreground):
            raise ValueError("Duplicate CTgene IDs")
        universe = set(row["gene_id"] for row in external)
        missing = [row for row in ct if row["gene_id"] not in universe]
        counts = Counter(row["gene_id"] for row in external)
        duplicates = [row for row in external if counts[row["gene_id"]] > 1]
        write_tsv(missing, output/f"{phase}_ct_missing_external.tsv", list(ct[0]))
        write_tsv(duplicates, output/f"{phase}_duplicate_rows.tsv", list(external[0]))

        # 2. Compute all counts from this same phase-specific measured universe.
        sets, universe_size, ct_count = prepare_membership(
            pathways, universe, foreground, config["min_size"], config["max_size"])
        if not sets:
            raise ValueError("No pathways satisfy the size bounds")
        metadata = [{key: value for key, value in p.items() if key != "genes"} for p in sets]
        p_values = hypergeometric(sets, universe_size, ct_count)
        adjusted = bh_adjust(p_values)  # Adjust BEFORE selecting pathways to display.
        results = [dict(row, hypergeometric_p=float(p), hypergeometric_bh=float(q))
                   for row, p, q in zip(metadata, p_values, adjusted)]
        write_tsv(results, output/f"{phase}_hypergeometric.tsv")
        summaries.append({"phase": phase, "ct_input": len(ct), "external_rows": len(external),
                          "universe_size": universe_size, "ct_in_universe": ct_count,
                          "missing_ct": len(missing), "duplicate_extra_rows": len(external)-universe_size,
                          "tested_pathways": len(sets), "pathways_with_ct": sum(p["ct_count"] > 0 for p in sets)})

        # 3. A score method runs only when its source column and meaning are set.
        column, kind = config.get("score_column"), config["score_kind"]
        score_field = config.get("score_field", "").format(phase=phase)
        if (column is None and not score_field) or kind is None:
            print(f"{phase}: membership analysis complete; score choice is unset.", flush=True)
            continue
        if not score_field and column not in range(2, 7):
            raise ValueError("score_column must be an original column number 2 through 6")
        field = score_field or f"source_col_{column}"
        values = transform_score([float(r[field]) for r in external], kind)
        scores = resolve_duplicate_scores([r["gene_id"] for r in external], values, config["duplicate_policy"])
        write_tsv([{"gene_id": gene, "analysis_score": value} for gene, value in scores.items()],
                  output/f"{phase}_analysis_scores.tsv")

        # 4. Each method compares the observed set against random sets of its size.
        for method in config["methods"]:
            print(f"{phase}: {method}, {config['permutations']} random sets per distinct size"
                  + (" (experimental)" if method == "maxmean" else ""), flush=True)
            test = permutation_test(scores, [p["genes"] for p in sets], method,
                                    config["permutations"], config["seed"] + (phase == "middle"), verbose=True)
            write_tsv([dict(row, **statistic) for row, statistic in zip(metadata, test)],
                      output/f"{phase}_{method}.tsv")
    write_tsv(summaries, output/"input_summary.tsv")
    print(f"Outputs: {output.resolve()}", flush=True)


if __name__ == "__main__":
    main()
