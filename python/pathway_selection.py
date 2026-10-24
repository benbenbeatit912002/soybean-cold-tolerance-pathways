"""Apply explicit Excel-style conditions to complete candidate method tables.

This module has no expected pathway IDs or target result counts. It supports
nominal or BH probabilities and separate biological rules for each phase.
"""
import math

METHODS = ("hypergeometric", "sumstat", "gsea")


def select_candidates(rows, config):
    field, alpha = config["probability_field"], config["alpha"]
    if field not in ("p_value", "p_bh") or not 0 < alpha <= 1:
        raise ValueError("Use a named probability field and 0 < alpha <= 1")
    if config["method_rule"] not in ("any", "all") or config["count_method"] not in METHODS:
        raise ValueError("Invalid method rule or count source")
    groups = {}
    for row in rows:
        phase, go, method = row["phase"], row["go_id"], row["method"]
        if phase not in ("short", "middle") or method not in METHODS or not go.strip():
            raise ValueError("Invalid phase, method or GO ID")
        methods = groups.setdefault((phase, go), {})
        if method in methods:
            raise ValueError("Duplicate phase/GO/method row")
        methods[method] = row
    audit = []
    for phase in ("short", "middle"):
        rule = config["phases"][phase]
        min_ct, min_ratio = rule["min_ct"], rule["min_ratio"]
        if not math.isfinite(min_ct) or min_ct < 0 or min_ct != int(min_ct) or not 0 <= min_ratio <= 1:
            raise ValueError("Invalid biological threshold")
        if rule["biological_rule"] not in ("count", "ratio", "and", "or"):
            raise ValueError("Invalid biological rule")
        for key in sorted(k for k in groups if k[0] == phase):
            methods = groups[key]
            if set(methods) != set(METHODS):
                raise ValueError(f"Incomplete method family for {key}")
            counts, sizes, values = {}, {}, {}
            for method, row in methods.items():
                n, ct, p = float(row["annotation_size"]), float(row["ct_count"]), float(row[field])
                if not all(map(math.isfinite, (n, ct, p))) or n < 1 or n != int(n) or ct < 0 or ct > n or ct != int(ct) or not 0 <= p <= 1:
                    raise ValueError(f"Invalid candidate metadata or probability for {key}")
                sizes[method], counts[method], values[method] = int(n), int(ct), p
            if len(set(sizes.values())) != 1:
                raise ValueError(f"Annotation sizes differ across methods for {key}")
            source = config["count_method"]
            count, size = counts[source], sizes[source]
            ratio = count / size
            count_pass, ratio_pass = count >= min_ct, ratio >= min_ratio
            biological = {"count": count_pass, "ratio": ratio_pass,
                          "and": count_pass and ratio_pass, "or": count_pass or ratio_pass}[rule["biological_rule"]]
            passing = [m for m in METHODS if values[m] < alpha]
            statistic = len(passing) >= (1 if config["method_rule"] == "any" else len(METHODS))
            audit.append(dict(phase=phase, go_id=key[1], description=methods[source]["description"],
                annotation_size=size, ct_count=count, ct_ratio=ratio, count_method=source,
                method_ct_counts_agree=int(len(set(counts.values())) == 1),
                probability_field=field, **{m+"_probability": values[m] for m in METHODS},
                count_pass=int(count_pass), ratio_pass=int(ratio_pass), biological_pass=int(biological),
                methods_passing=";".join(passing), statistic_pass=int(statistic),
                selected=int(biological and statistic), selection_basis=config["label"]))
    if not audit or {r["phase"] for r in audit} != {"short", "middle"}:
        raise ValueError("Candidate input must contain both phases")
    return audit
