"""Recalculate short scores from expression, correspondence, lengths and SNP blocks.

Usage: python scripts/reproduce_historical_short.py ORIGINAL_PROJECT OUTPUT_DIRECTORY
No saved gene_p_len, corrected probabilities, or reference scores are read.
Exit 2 means the reconstructed historical procedure left nonfinite scores;
outputs and diagnostics are still written. Verification is a separate command.
"""
import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
from ctpathway import write_tsv
from historical_short import (reconstruct_input, historical_snp_imputation, fit_historical_short,
                              historical_score_export, repair_zero_hypothesis)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--zero-repair-hypothesis", action="store_true",
                        help="test whole-genome bin mean, last-phase lambda, and 7-digit repair; provenance unconfirmed")
    args = parser.parse_args()
    rows, exclusions, metadata = reconstruct_input(args.source_root)
    counts, audit = historical_snp_imputation(rows)
    for row, count in zip(rows, counts):
        row["historical_snp"] = int(count)
    exponent, adjusted, scores, grid = fit_historical_short([r["raw_short"] for r in rows], counts)
    args.output.mkdir(parents=True, exist_ok=True)
    if args.zero_repair_hypothesis:
        counts, adjusted, repair_audit, middle_grid = repair_zero_hypothesis(rows, counts, adjusted, args.source_root)
        if repair_audit:
            write_tsv(repair_audit, args.output/"zero_repair_hypothesis.tsv")
        write_tsv(middle_grid, args.output/"middle_lambda_grid_for_hypothesis.tsv")
    unresolved = []
    for row, count, probability in zip(rows, counts, adjusted):
        saved_precision, score, display = historical_score_export(probability)
        row.update(snp=int(count), adjusted_short=float(probability),
                   adjusted_short_saved_precision=saved_precision,
                   score_short_unrounded=score, score_short=display,
                   calculation_status="finite" if math.isfinite(score) else "unresolved_zero_probability")
        if not math.isfinite(score):
            unresolved.append(dict(row))
    write_tsv(rows, args.output/"recalculated_short.tsv")
    write_tsv(exclusions, args.output/"mapping_exclusions.tsv",
              ["expression_row", "original_gene_id", "reason"])
    write_tsv(audit, args.output/"historical_snp_imputation.tsv",
              ["source_row", "gene_id", "bin_start", "bin_end_exclusive", "group_records",
               "mean_including_zeros", "snp_before", "snp_after"])
    write_tsv(grid, args.output/"short_lambda_grid.tsv")
    write_tsv(unresolved, args.output/"unresolved_rows.tsv", list(rows[0]))
    metadata.update(lambda_short=exponent, zero_markers_before=sum(r["snp_markers"] == 0 for r in rows),
                    zero_snp_after=int(sum(counts == 0)), nonfinite_scores=len(unresolved),
                    arithmetic="1-(1-p)^(SNP^lambda); ordinary floating-point subtraction",
                    zero_policy="original_bin_mean_including_zeros_no_fallback",
                    export_rule="15-significant-digit probability; zero score to 1e-16; inferred saved decimal format",
                    zero_repair_hypothesis=args.zero_repair_hypothesis,
                    repair_provenance_confirmed=False, reference_scores_read=False, python=sys.version)
    (args.output/"recalculation.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Recalculated {len(rows)} records; lambda={exponent}; unresolved={len(unresolved)}")
    if args.zero_repair_hypothesis:
        print("HYPOTHESIS: matching saved values does not confirm the historical repair procedure.")
    return 2 if unresolved else 0


if __name__ == "__main__":
    sys.exit(main())
