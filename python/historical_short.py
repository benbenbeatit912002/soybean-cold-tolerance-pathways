"""Reconstruct the archived short-term calculation from four upstream files.

This module never reads corrected probabilities or published/reference scores.
It deliberately follows the archived arithmetic and zero-SNP rule, including
their limitations, so that discrepancies remain visible to the verifier.
"""
from collections import defaultdict
import csv
import hashlib
import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from ctpathway import canonical_id

SOURCE_FILES = {
    "expression": "研究法/大豆基因/數據/外部資料匯入/耐寒資料庫0 1 24hr(轉換後)培修版.txt",
    "correspondence": "研究法/大豆基因/大豆基因 外部/Glyma_genes_correspondence_final.txt",
    "gene_length": "研究法/大豆基因/數據/外部基因校正/Glyma2.0.txt",
    "snp_markers": "研究法/大豆基因/數據/外部基因校正/snps_allgenes_20kb.set.txt",
}


def resolve_source_files(source_root):
    """Read a portable manifest, or retain support for the original archive.

    A manifest is authoritative: missing or changed packaged files are errors,
    never a reason to fall back to a different copy on the author's computer.
    """
    root = Path(source_root)
    manifest = root / "inputs.tsv"
    if not manifest.is_file():
        paths = {name: root / relative for name, relative in SOURCE_FILES.items()}
    else:
        with manifest.open(encoding="utf-8-sig", newline="") as handle:
            entries = list(csv.DictReader(handle, delimiter="\t"))
        if len(entries) != len(SOURCE_FILES) or {r["source"] for r in entries} != set(SOURCE_FILES):
            raise ValueError("inputs.tsv must contain each of the four upstream sources once")
        paths = {}
        for entry in entries:
            path = root / entry["file"]
            contents = path.read_bytes()
            if (hashlib.sha256(contents).hexdigest() != entry["sha256"] or
                    hashlib.md5(contents).hexdigest() != entry["md5"]):
                raise ValueError(f"Upstream checksum mismatch: {path}")
            paths[entry["source"]] = path
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing upstream inputs: " + ", ".join(missing))
    return paths


def read_table(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))[1:]


def count_snp_markers(path):
    """Count AX markers between each gene header and END; retain empty blocks."""
    counts = {}
    gene = None
    count = 0
    with Path(path).open(encoding="utf-8-sig") as handle:
        for line in handle:
            value = line.strip()
            if not value:
                continue
            if value.startswith("G"):
                if gene is not None or value in counts:
                    raise ValueError("Missing END or duplicate SNP gene block")
                gene, count = value, 0
            elif value.startswith("A") and gene is not None:
                count += 1
            elif value == "END" and gene is not None:
                counts[gene] = count
                gene = None
            else:
                raise ValueError(f"Unexpected SNP block entry: {value}")
    if gene is not None:
        raise ValueError("Unclosed SNP block")
    return counts


def reconstruct_input(source_root):
    """Match unique correspondence ROWS, exactly as the original R loop did."""
    sources = resolve_source_files(source_root)
    correspondence = read_table(sources["correspondence"])
    by_alias = defaultdict(set)
    for index, row in enumerate(correspondence):
        for alias in row[:2]:
            by_alias[alias].add(index)
    lengths = {}
    for row in read_table(sources["gene_length"]):
        if row[0] in lengths:
            raise ValueError("Duplicate gene-length key")
        lengths[row[0]] = int(row[4])
    snp_counts = count_snp_markers(sources["snp_markers"])
    rows = []
    excluded = []
    expression = read_table(sources["expression"])
    for line, source in enumerate(expression, start=2):
        matches = by_alias[source[2]]
        reason = ""
        if not matches:
            reason = "no_correspondence"
        elif len(matches) != 1:
            reason = "ambiguous_correspondence_rows"
        else:
            gene = correspondence[next(iter(matches))][2]
            if not gene:
                reason = "blank_target"
        if reason:
            excluded.append({"expression_row": line, "original_gene_id": source[2], "reason": reason})
            continue
        # Columns 5 and 8 are supplied expression p-values, not fold changes.
        rows.append({"source_row": len(rows)+2, "expression_row": line,
                     "original_gene_id": source[2], "gene_id": canonical_id(gene),
                     "raw_short": source[4], "raw_middle": source[7],
                     "gene_length": lengths[gene], "snp_markers": snp_counts[gene]})
    metadata = {"expression_records": len(expression), "mapped_records": len(rows),
                "unique_genes": len({r["gene_id"] for r in rows}),
                "source_sha256": {name: hashlib.sha256(p.read_bytes()).hexdigest() for name, p in sources.items()},
                "source_paths": {name: str(p) for name, p in sources.items()}}
    return rows, excluded, metadata


def historical_snp_imputation(rows):
    """Replicate old 2,000-bp bins: the mean INCLUDES zeros, with no fallback.

    The original seq(1,100000,by=2000) ends at 98001. Its last two loop
    iterations have no usable upper boundary. They therefore change no records.
    """
    values = np.array([r["snp_markers"] for r in rows], dtype=float)
    lengths = np.array([r["gene_length"] for r in rows], dtype=float)
    audit = []
    for start in range(1, 98001, 2000):
        members = (lengths >= start) & (lengths < start+2000)
        zeros = np.flatnonzero(members & (values == 0))
        if not zeros.size:
            continue
        mean = float(np.mean(values[members]))
        replacement = round(mean)  # R round(): nearest integer, ties to even.
        for index in zeros:
            audit.append({"source_row": rows[index]["source_row"], "gene_id": rows[index]["gene_id"],
                          "bin_start": start, "bin_end_exclusive": start+2000,
                          "group_records": int(members.sum()), "mean_including_zeros": mean,
                          "snp_before": 0, "snp_after": replacement})
            values[index] = replacement
    return values, audit


def historical_fosco(probabilities, snp, exponent):
    """Keep ordinary subtraction: stable log1p/expm1 changes archived scores."""
    return np.array([1-math.pow(1-float(p), math.pow(float(n), exponent))
                     for p, n in zip(probabilities, snp)])


def fit_historical_short(probabilities, snp):
    p = np.asarray(probabilities, dtype=float)
    n = np.asarray(snp, dtype=float)
    if p.shape != n.shape or p.ndim != 1 or len(p) < 2:
        raise ValueError("Expected matching nontrivial probability and SNP vectors")
    if not np.isfinite(p).all() or not ((p > 0) & (p <= 1)).all():
        raise ValueError("Raw probabilities must be finite and in (0,1]")
    if not np.isfinite(n).all() or (n < 0).any() or (n != np.rint(n)).any():
        raise ValueError("SNP counts must be finite nonnegative integers")
    grid = []
    for step in range(1, 101):
        exponent = step/100
        adjusted = historical_fosco(p, n, exponent)
        rho = float(spearmanr(n, adjusted).statistic)
        if not math.isfinite(rho):
            raise ValueError("Cannot fit lambda when the Spearman correlation is undefined")
        grid.append({"lambda": exponent, "absolute_spearman": abs(rho)})
    best = min(grid, key=lambda row: row["absolute_spearman"])["lambda"]
    adjusted = historical_fosco(p, n, best)
    scores = np.array([-math.log10(value) if value > 0 else math.inf for value in adjusted])
    return best, adjusted, scores, grid


def historical_score_export(adjusted):
    """Recreate the precision pattern observed in every finite saved short score.

    R write.table exported corrected probabilities to 15 significant digits.
    The saved log file replaces zero scores by 1e-16 and uses a mixed fixed /
    scientific decimal format. The export rules are inferred from the files;
    the original spreadsheet/export command has not been recovered.
    """
    probability = float(format(adjusted, ".15g"))
    score = -math.log10(probability) if probability > 0 else math.inf
    if score == 0:
        score = 1e-16
    if not math.isfinite(score):
        display = "Inf"
    elif score < 1e-4:
        display = format(score, ".2E")
    else:
        display = format(score, ".8f" if score >= 10 else ".9f").rstrip("0").rstrip(".")
    return probability, score, display


def repair_zero_hypothesis(rows, counts, adjusted, source_root):
    """Test an unconfirmed historical repair, without looking up saved scores.

    For zero counts left by the old rule, use positive SNP donors from the same
    bin in the whole genome. Recalculate with the *last phase's* fitted lambda
    (middle) and 7 significant digits, consistent with an R console value.
    This is an explicit reconstruction hypothesis, not an established method.
    """
    repaired_counts = counts.copy()
    repaired_adjusted = adjusted.copy()
    last_lambda, _, _, middle_grid = fit_historical_short([r["raw_middle"] for r in rows], counts)
    paths = resolve_source_files(source_root)
    genome_counts = count_snp_markers(paths["snp_markers"])
    genome = [(r[0], int(r[4]), genome_counts.get(r[0])) for r in read_table(paths["gene_length"])]
    audit = []
    for index in np.flatnonzero(counts == 0):
        row = rows[index]
        start = 1+2000*((row["gene_length"]-1)//2000)
        donors = [(gene, n) for gene, length, n in genome
                  if start <= length < start+2000 and n is not None and n > 0]
        if not donors:
            continue
        mean = sum(n for _, n in donors)/len(donors)
        replacement = round(mean)
        probability = historical_fosco([row["raw_short"]], [replacement], last_lambda)[0]
        rounded = float(format(probability, ".7g"))
        repaired_counts[index] = replacement
        repaired_adjusted[index] = rounded
        audit.append({"source_row": row["source_row"], "gene_id": row["gene_id"],
                      "bin_start": start, "bin_end_exclusive": start+2000,
                      "donor_genes": ";".join(gene for gene, _ in donors),
                      "donor_counts": ";".join(str(n) for _, n in donors),
                      "donor_mean": mean, "snp_after": replacement,
                      "lambda_from_middle_grid": last_lambda,
                      "adjusted_before_rounding": float(probability), "adjusted_7_significant_digits": rounded,
                      "evidence_status": "unconfirmed_reconstruction_hypothesis"})
    return repaired_counts, repaired_adjusted, audit, middle_grid
