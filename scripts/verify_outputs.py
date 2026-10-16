"""Verify staged hashes and completed analysis outputs; emit a compact audit."""
import csv
import hashlib
import json
import math
from pathlib import Path


def table(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "reports/source_manifest.json").read_text(encoding="utf-8"))
    for item in manifest:
        assert hashlib.sha256((root/item["staged"]).read_bytes()).hexdigest() == item["staged_sha256"], item["staged"]
        source = root.parent / item["source"]
        if source.exists():
            assert hashlib.sha256(source.read_bytes()).hexdigest() == item["source_sha256"], item["source"]
    summary = {"verified_staged_hashes": len(manifest), "runs": {}}
    for name in ("external_membership", "score5_diagnostic", "score6_diagnostic"):
        folder = root / "results" / name
        info = {"qc": table(folder / "input_summary.tsv"), "tests": []}
        for phase in ("short", "middle"):
            for method in (("hypergeometric",) if name == "external_membership" else ("sumstat", "gsea")):
                rows = table(folder/f"{phase}_{method}.tsv")
                pcol, qcol = ("hypergeometric_p", "hypergeometric_bh") if method == "hypergeometric" else ("p_value", "p_bh")
                pv = [float(r[pcol]) for r in rows]
                qv = [float(r[qcol]) for r in rows]
                assert rows and len({r["go_id"] for r in rows}) == len(rows)
                assert all(math.isfinite(p) and 0 <= p <= q <= 1 for p, q in zip(pv, qv))
                if method != "hypergeometric":
                    assert min(pv) >= 1/200
                    assert all(math.isfinite(float(r["statistic"])) for r in rows)
                info["tests"].append({"phase": phase, "method": method, "pathways": len(rows),
                                      "min_p": min(pv), "min_bh": min(qv), "bh_below_0_05": sum(q < .05 for q in qv)})
        summary["runs"][name] = info
    (root/"reports/output_verification.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
