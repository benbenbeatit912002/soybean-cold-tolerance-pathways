"""Compare fresh scalar/vector calculations with saved short-term validation data."""
import csv
import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results/short_derivation_2026-10-21"
def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

raw = read(ROOT/"data/reference/historical_input.tsv")
expected = read(ROOT/"data/reference/paper_short.tsv")
r_output = read(OUT/"r_naive_short.tsv")
p = np.array([float(row["raw_short"]) for row in raw])
snp = np.array([float(row["snp"]) for row in raw])
legacy = np.array([float(row["reference_adjusted_short"]) for row in expected])
score = np.array([float(row["reference_score_short"]) for row in expected])
candidates = {
    "r_naive": np.array([float(row["adjusted_short"]) for row in r_output]),
    "python_scalar_pow": np.array([1-math.pow(1-x, math.pow(n, .01)) for x, n in zip(p, snp)]),
    "numpy_power": 1-np.power(1-p, np.power(snp, .01)),
    "stable": -np.expm1(np.power(snp, .01)*np.log1p(-p)),
    "saved_corrected_probability_log_only": legacy,
}
summary = {}
rows = []
for name, probability in candidates.items():
    with np.errstate(divide="ignore"):
        calculated = -np.log10(probability)
    delta = calculated-score
    positive = snp > 0
    summary[name] = {
        "probability_max_absolute_error_positive_snp": float(np.max(np.abs(probability[positive]-legacy[positive]))),
        "score_max_absolute_error_positive_snp": float(np.max(np.abs(delta[positive]))),
        "score_rows_within_5e_8": int(np.sum(np.abs(delta) <= 5e-8)),
        "score_rows_outside_5e_8_positive_snp": int(np.sum(np.abs(delta[positive]) > 5e-8)),
        "probability_bitwise_differences_from_r": int(np.sum(probability != candidates["r_naive"])),
    }
    for i in np.flatnonzero(np.abs(delta) > 5e-8):
        rows.append({"method": name, "source_row": raw[i]["source_row"], "gene_id": raw[i]["gene_id"],
                     "raw_short": p[i], "snp": snp[i], "calculated_adjusted": probability[i],
                     "saved_adjusted": legacy[i], "calculated_score": calculated[i],
                     "saved_score": score[i], "difference": delta[i]})
with (OUT/"numeric_differences.tsv").open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
(OUT/"numeric_summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
print(json.dumps(summary, indent=2))
