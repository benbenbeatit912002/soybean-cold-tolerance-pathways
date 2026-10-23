"""Readable statistical functions for soybean cold-tolerance pathways.

The matching R implementation is R/core.R. Inputs are already mapped to one
genome version by scripts/prepare_inputs.py. A gene's prior CT score is never
added to its external score. See docs/METHODS.md for the statistical assumptions.
"""
from collections import defaultdict
import csv
from pathlib import Path
import re

import numpy as np
from scipy import stats


def canonical_id(value):
    """Normalize g/G spelling; this is NOT a genome-version conversion."""
    return re.sub(r"^Glyma\.(\d{2})[gG]", r"Glyma.\1G", str(value).strip())


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(rows, path, columns=None):
    """Always write column names, including when the result has no rows."""
    if columns is None:
        if not rows:
            raise ValueError("Supply columns for an empty table")
        columns = list(rows[0])
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_gmt(path):
    """A GMT row is: GO ID, description, gene 1, gene 2, ... ."""
    pathways = []
    with Path(path).open(encoding="utf-8-sig") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3:
                raise ValueError("Invalid GMT row")
            pathways.append({"go_id": fields[0], "description": fields[1],
                             "genes": list(dict.fromkeys(canonical_id(x) for x in fields[2:]))})
    if not pathways or len({p["go_id"] for p in pathways}) != len(pathways):
        raise ValueError("Empty GMT or duplicate pathway IDs")
    return pathways


def validate_values(values, kind):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not values.size or not np.isfinite(values).all():
        raise ValueError("Scores must be a nonempty finite numeric vector")
    if kind == "pvalue" and np.any((values <= 0) | (values > 1)):
        raise ValueError("pvalue requires 0 < value <= 1; resolve zeros upstream")
    if kind == "fold_change" and np.any(values <= 0):
        raise ValueError("fold_change requires positive ratios")
    if kind == "score" and np.any(values < 0):
        raise ValueError("score requires nonnegative, larger-is-stronger values")
    if kind not in ("pvalue", "fold_change", "score"):
        raise ValueError("Unknown score kind")
    return values


def transform_score(values, kind):
    """Return a larger-is-stronger score; leave supplied scores unchanged."""
    values = validate_values(values, kind)
    if kind == "pvalue":
        return -np.log10(values)
    if kind == "fold_change":
        return np.abs(np.log(values))
    return values.copy()


def resolve_duplicate_scores(ids, values, policy="error"):
    """Retain one score per gene. Averaging must be explicitly requested."""
    if policy not in ("error", "mean") or len(ids) != len(values):
        raise ValueError("Invalid duplicate policy or mismatched input lengths")
    groups = defaultdict(list)
    for gene, value in zip(ids, values):
        gene = canonical_id(gene)
        if not gene:
            raise ValueError("Missing gene ID")
        groups[gene].append(value)
    repeated = any(len(group) > 1 for group in groups.values())
    if repeated and policy == "error":
        raise ValueError("Duplicate gene IDs: resolve provenance or explicitly choose mean")
    # R's mean policy sorts gene names. Use the same ordering here.
    order = sorted(groups) if repeated else groups.keys()
    return {gene: float(np.mean(groups[gene])) for gene in order}


def prepare_membership(pathways, universe, foreground, min_size=6, max_size=1500):
    """Restrict every pathway and the CT list to the same measured universe."""
    universe = {canonical_id(g) for g in universe}
    foreground = {canonical_id(g) for g in foreground} & universe
    prepared = []
    for pathway in pathways:
        genes = [g for g in pathway["genes"] if g in universe]
        if not min_size <= len(genes) <= max_size:
            continue
        ct_count = len(set(genes) & foreground)
        prepared.append({"go_id": pathway["go_id"], "description": pathway["description"],
                         "annotation_size": len(pathway["genes"]), "tested_size": len(genes),
                         "ct_count": ct_count, "ct_ratio": ct_count/len(genes), "genes": genes})
    return prepared, len(universe), len(foreground)


def hypergeometric(pathways, universe_size, ct_in_universe):
    """P(overlap >= observed), computed directly as an upper tail."""
    return np.asarray([stats.hypergeom.sf(p["ct_count"] - 1, universe_size,
                                         p["tested_size"], ct_in_universe) for p in pathways])


def bh_adjust(p_values):
    """Benjamini-Hochberg adjustment across the entire supplied test family."""
    p_values = np.asarray(p_values, dtype=float)
    if not p_values.size:
        return p_values.copy()
    if np.any(~np.isfinite(p_values) | (p_values < 0) | (p_values > 1)):
        raise ValueError("Invalid p-values")
    order = np.argsort(p_values, kind="stable")
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty_like(p_values)
    result[order] = np.minimum(adjusted, 1)
    return result


def monte_carlo_p(observed, null):
    """Include ties and the observed result: (exceedances + 1)/(draws + 1)."""
    null = np.asarray(null, dtype=float)
    if not np.isfinite(observed) or not null.size or not np.isfinite(null).all():
        raise ValueError("Invalid permutation statistics")
    tolerance = 1e-12 * max(1, abs(observed))
    return float((1 + np.count_nonzero(null >= observed-tolerance)) / (len(null)+1))


def maxmean_stat(z_scores):
    """Larger of the positive and negative contributions per SET member."""
    z_scores = np.asarray(z_scores, dtype=float)
    return float(max(np.maximum(z_scores, 0).mean(), np.maximum(-z_scores, 0).mean()))


def gsea_es(ranked_scores, hits):
    """Positive-tail weighted running sum. Python hit positions are zero-based.

    Only hit positions need evaluation: the running sum decreases between hits.
    This is a GSEA-style score, not Broad GSEA's NES or empirical FDR.
    """
    ranked_scores = np.asarray(ranked_scores, dtype=float)
    hits = np.unique(np.asarray(hits, dtype=int))
    n_genes, n_hits = len(ranked_scores), len(hits)
    if not 0 < n_hits < n_genes or hits[0] < 0 or hits[-1] >= n_genes:
        raise ValueError("GSEA needs 0 < set size < universe size and valid indices")
    weights = np.abs(ranked_scores[hits])
    if weights.sum() == 0:
        weights = np.ones(n_hits)
    cumulative_hits = np.cumsum(weights)/weights.sum()
    misses_before_hit = hits - np.arange(n_hits)
    running_sum = cumulative_hits - misses_before_hit/(n_genes-n_hits)
    return float(max(0, running_sum.max()))


def permutation_test(scores, members, method, permutations=10000, seed=20261015, verbose=False):
    """Random gene-set test, using one null distribution per distinct set size.

    Sampling is WITHOUT replacement within each draw, and restarts for the next
    draw. Python and R use different random-number generators, so their Monte
    Carlo estimates need not be identical even for the same numeric seed.
    """
    if not isinstance(permutations, int) or isinstance(permutations, bool) or permutations < 1:
        raise ValueError("permutations must be a positive integer")
    if method not in ("sumstat", "gsea", "maxmean"):
        raise ValueError("Unknown method")
    gene_ids = list(scores)
    values = validate_values(list(scores.values()), "score")
    if method == "gsea":
        gene_ids = sorted(gene_ids, key=lambda gene: (-scores[gene], gene))
        values = np.array([scores[gene] for gene in gene_ids])
    if method == "maxmean":
        deviation = values.std(ddof=1)  # R sd() also uses n-1.
        if not np.isfinite(deviation) or deviation == 0:
            raise ValueError("MaxMean requires variable scores")
        values = (values-values.mean())/deviation
    positions = {gene: i for i, gene in enumerate(gene_ids)}
    indices = [np.array([positions[gene] for gene in genes], dtype=int) for genes in members]
    if any(not len(index) or len(set(index)) != len(index) for index in indices):
        raise ValueError("Pathways need distinct measured genes")

    def statistic(index):
        if method == "sumstat":
            return float(values[index].sum())
        if method == "maxmean":
            return maxmean_stat(values[index])
        return gsea_es(values, index)

    observed = np.array([statistic(index) for index in indices])
    sizes = np.array([len(index) for index in indices])
    distinct_sizes = sorted(set(sizes))
    p_values = np.empty(len(indices))
    rng = np.random.default_rng(seed)
    for number, size in enumerate(distinct_sizes, start=1):
        null = np.empty(permutations)
        for draw in range(permutations):
            random_genes = rng.choice(len(values), size=int(size), replace=False)
            null[draw] = statistic(random_genes)
        for pathway_index in np.flatnonzero(sizes == size):
            p_values[pathway_index] = monte_carlo_p(observed[pathway_index], null)
        if verbose and (number % 50 == 0 or number == len(distinct_sizes)):
            print(f"{method}: completed {number}/{len(distinct_sizes)} distinct pathway sizes", flush=True)
    adjusted = bh_adjust(p_values)
    return [{"statistic": float(value), "p_value": float(p), "p_bh": float(q),
             "null_exceedances": int(round(p*(permutations+1)-1)),
             "permutations": permutations, "seed": seed}
            for value, p, q in zip(observed, p_values, adjusted)]


def match_short_reference(rows, reference):
    """Match a saved short-score reference by source row, checking input identity.

    A gene can have multiple records. Never join on gene ID alone or average
    these reference scores. Return original reference strings without rounding.
    """
    identity = {"source_row", "gene_id", "raw_short", "gene_length", "snp"}
    saved = {"reference_adjusted_short", "reference_score_short",
             "reference_snp_short", "reference_log_row"}
    if (not rows or not reference or any(not identity.issubset(row) for row in rows)
            or any(not (identity | saved).issubset(row) for row in reference)):
        raise ValueError("Short reference and input require source-row identity and saved score columns")

    def positive_integer(value):
        number = float(value)
        if not np.isfinite(number) or number < 1 or number != np.floor(number):
            raise ValueError("Reference row numbers and SNP counts must be positive integers")
        return int(number)

    keys = [positive_integer(row["source_row"]) for row in rows]
    reference_keys = [positive_integer(row["source_row"]) for row in reference]
    if (len(set(keys)) != len(keys) or len(set(reference_keys)) != len(reference_keys)
            or set(keys) != set(reference_keys)):
        raise ValueError("Short reference source rows must match input exactly, without duplicates")
    lookup = dict(zip(reference_keys, reference))
    matched = [lookup[key] for key in keys]
    for row, old in zip(rows, matched):
        gene = canonical_id(row["gene_id"])
        if not gene or gene != canonical_id(old["gene_id"]):
            raise ValueError("Short reference gene identity differs from input")
        for field in ("raw_short", "gene_length", "snp"):
            value, original = float(row[field]), float(old[field])
            if not np.isfinite(value) or value != original:
                raise ValueError(f"Short reference {field} differs from input")
        positive_integer(old["reference_snp_short"])
        positive_integer(old["reference_log_row"])
    validate_values([float(row["reference_adjusted_short"]) for row in matched], "pvalue")
    validate_values([float(row["reference_score_short"]) for row in matched], "score")
    return matched


def impute_snp_counts(snp_counts, gene_lengths=None, bin_width=2000):
    """Fill zeros from positive input counts in the same length bin, else globally.

    Bins are [1, 2001), [2001, 4001), ... with no upper length limit.
    Donors are frozen before filling: an imputed count never becomes a donor.
    Round to the nearest integer, with ties to even, matching R's round().
    """
    counts = np.asarray(snp_counts, dtype=float)
    if (counts.ndim != 1 or not counts.size or
            np.any(~np.isfinite(counts) | (counts < 0) | (counts != np.floor(counts)))):
        raise ValueError("SNP counts must be finite nonnegative integers")
    if not np.isfinite(bin_width) or bin_width < 1 or bin_width != int(bin_width):
        raise ValueError("bin_width must be a positive integer")
    positive = counts > 0
    if not positive.any():
        raise ValueError("Cannot impute SNP counts: no positive input donors")
    if gene_lengths is None:
        if not positive.all():
            raise ValueError("gene_length is required when SNP counts contain zero")
        return {"values": counts.copy(), "audit": []}
    lengths = np.asarray(gene_lengths, dtype=float)
    if (lengths.shape != counts.shape or
            np.any(~np.isfinite(lengths) | (lengths < 1) | (lengths != np.floor(lengths)))):
        raise ValueError("Gene lengths must be finite positive integers matching SNP counts")
    bins = np.floor((lengths - 1) / bin_width)
    filled = counts.copy()
    audit = []
    for group in np.unique(bins[~positive]):
        donors = counts[positive & (bins == group)]
        source = "group_nonzero_mean"
        if not donors.size:
            donors = counts[positive]
            source = "global_nonzero_mean"
        donor_mean = float(np.mean(donors))
        replacement = int(np.rint(donor_mean))
        for index in np.flatnonzero((bins == group) & ~positive):
            filled[index] = replacement
            audit.append({"row_number": int(index + 1),
                          "bin_start": int(group * bin_width + 1),
                          "bin_end_exclusive": int((group + 1) * bin_width + 1),
                          "original_snp": 0, "imputed_snp": replacement,
                          "donor_count": int(donors.size), "donor_mean": donor_mean,
                          "source": source})
    return {"values": filled, "audit": sorted(audit, key=lambda row: row["row_number"])}


def fosco(probabilities, snp_counts, tuning_parameter):
    """Stable 1-(1-p)^(M^lambda); requires confirmed probabilities or pseudo-scores."""
    probabilities = validate_values(probabilities, "pvalue")
    snp_counts = np.asarray(snp_counts, dtype=float)
    if snp_counts.shape != probabilities.shape or np.any(~np.isfinite(snp_counts) | (snp_counts <= 0)):
        raise ValueError("FOSCO needs a positive SNP count for every score")
    if not np.isfinite(tuning_parameter) or tuning_parameter < 0:
        raise ValueError("Invalid lambda")
    with np.errstate(divide="ignore"):
        return -np.expm1(snp_counts**tuning_parameter * np.log1p(-probabilities))


def fit_fosco(probabilities, snp_counts, grid=None):
    """Select the smallest grid lambda minimizing absolute Spearman correlation."""
    grid = np.arange(1, 101)/100 if grid is None else np.asarray(grid, dtype=float)
    if not grid.size or np.any(~np.isfinite(grid) | (grid < 0)):
        raise ValueError("Invalid lambda grid")
    objective = []
    for tuning_parameter in grid:
        adjusted = fosco(probabilities, snp_counts, tuning_parameter)
        if len(set(snp_counts)) < 2 or len(set(adjusted)) < 2:
            objective.append(float("inf"))
        else:
            objective.append(abs(float(stats.spearmanr(snp_counts, adjusted).statistic)))
    if not np.isfinite(objective).any():
        raise ValueError("Cannot fit FOSCO: constant scores or SNP counts")
    selected = int(np.argmin(objective))
    return {"lambda": float(grid[selected]), "adjusted": fosco(probabilities, snp_counts, grid[selected]),
            "grid": [{"lambda": float(lam), "absolute_spearman": float(value)}
                     for lam, value in zip(grid, objective)]}
