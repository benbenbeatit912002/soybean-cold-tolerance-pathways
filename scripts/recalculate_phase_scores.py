"""Recalculate both phases from upstream files for enrichment analysis.

The explicit hypotheses option enables the documented zero-SNP repair and
replaces a cancellation-induced infinite middle score by the calculated maximum
finite score rounded to seven significant digits. No reference scores are read.
"""
import argparse
import json
import math
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
from ctpathway import write_tsv
from historical_short import (reconstruct_input, historical_snp_imputation,
    fit_historical_short, historical_fosco, historical_score_export, repair_zero_hypothesis)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--use-reconstruction-hypotheses", action="store_true")
    args = parser.parse_args()
    rows, excluded, metadata = reconstruct_input(args.source_root)
    counts, imputation = historical_snp_imputation(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    for row, count in zip(rows, counts):
        row["historical_snp"] = int(count)
    fits = {phase: fit_historical_short([r[f"raw_{phase}"] for r in rows], counts)
            for phase in ("short", "middle")}
    final_counts = counts.copy()
    repair_audit = []
    if args.use_reconstruction_hypotheses:
        final_counts, short_adjusted, repair_audit, _ = repair_zero_hypothesis(
            rows, counts, fits["short"][1], args.source_root)
    else:
        short_adjusted = fits["short"][1]
    exceptions = []
    for phase in ("short", "middle"):
        exponent, adjusted, _, grid = fits[phase]
        adjusted = short_adjusted.copy() if phase == "short" else adjusted.copy()
        if phase == "middle" and args.use_reconstruction_hypotheses:
            for i, (old, new) in enumerate(zip(counts, final_counts)):
                if old == 0 and new > 0:
                    adjusted[i] = float(format(historical_fosco([rows[i]["raw_middle"]], [new], exponent)[0], ".7g"))
        exported = [historical_score_export(value) for value in adjusted]
        finite_max = max(score for _, score, _ in exported if math.isfinite(score))
        for i, (row, value, (serialized, score, display)) in enumerate(zip(rows, adjusted, exported)):
            status = "recalculated"
            if counts[i] == 0 and final_counts[i] > 0:
                status = "zero_snp_repair_hypothesis"
            if not math.isfinite(score) and final_counts[i] > 0 and args.use_reconstruction_hypotheses:
                # The archived finite replacement is consistent with printing
                # the largest finite score using R's default seven digits.
                display = format(finite_max, ".7g")
                score = float(display)
                status = "maximum_finite_score_repair_hypothesis"
            if not math.isfinite(score):
                status = "unresolved_nonfinite"
            row.update(snp=int(final_counts[i]))
            row[f"adjusted_{phase}"] = float(value)
            row[f"adjusted_{phase}_serialized"] = serialized
            row[f"score_{phase}_unrounded"] = score
            row[f"score_{phase}"] = display
            row[f"status_{phase}"] = status
            if status != "recalculated":
                exceptions.append({"phase": phase, "source_row": row["source_row"], "gene_id": row["gene_id"],
                    "status": status, "historical_snp": int(counts[i]), "final_snp": int(final_counts[i]),
                    "adjusted_probability": float(value), "score": display, "calculated_finite_max": finite_max})
        write_tsv(grid, args.output/f"{phase}_lambda_grid.tsv")
    write_tsv(rows, args.output/"recalculated_scores.tsv")
    write_tsv(exceptions, args.output/"exception_audit.tsv")
    if repair_audit:
        write_tsv(repair_audit, args.output/"snp_repair_donors.tsv")
    write_tsv(imputation, args.output/"historical_snp_imputation.tsv")
    metadata.update(hypotheses_enabled=args.use_reconstruction_hypotheses,
                    historical_repairs_confirmed=False, reference_scores_read=False,
                    lambda_short=fits["short"][0], lambda_middle=fits["middle"][0])
    (args.output/"score_manifest.json").write_text(json.dumps(metadata, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(f"Both phases recalculated: {len(rows)} records; lambdas {fits['short'][0]}, {fits['middle'][0]}")
    if any(r[f"status_{phase}"] == "unresolved_nonfinite" for r in rows for phase in fits):
        return 2
    return 0

if __name__ == "__main__":
    sys.exit(main())
