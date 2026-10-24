"""Known-set checks for overlap, selection provenance, and invalid mixed inputs."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from presoynet import METHODS, jaccard_edges, select_pathways, filter_edges


class CrosstalkTests(unittest.TestCase):
    def test_historical_rounding_boundary(self):
        # Five shared genes, 1,122 genes in the union: rounding changes inclusion.
        members = {"a": set(range(50)), "b": set(range(45,1122))}
        edges, cutoff = filter_edges(members, .8, "historical-code", "short")
        self.assertLess(edges[0]["jaccard"], cutoff)
        self.assertEqual(edges[0]["filter_jaccard"], cutoff)
        self.assertEqual(edges[0]["retained"], 1)
        edges, _ = filter_edges(members, .5, "historical-code", "middle")
        self.assertEqual(edges[0]["retained"], 0)

    def test_overlap_ties_and_isolates(self):
        members = {"a": {"x", "y"}, "b": {"y", "z"}, "c": {"x", "z"}, "d": {"q"}}
        edges, cutoff = jaccard_edges(members, .5)
        self.assertEqual(len(edges), 6)
        self.assertEqual(cutoff, 1/3)
        self.assertEqual(sum(e["retained"] for e in edges), 3)  # Keep all three cutoff ties.
        self.assertTrue(all(e["source"] < e["target"] for e in edges))
        self.assertTrue(all(not e["retained"] for e in edges if e["target"] == "d"))
        self.assertEqual(sum(e["retained"] for e in jaccard_edges(members, 0)[0]), 0)
        self.assertEqual(sum(e["retained"] for e in jaccard_edges(members, 1)[0]), 3)

    def test_empty_single_and_invalid(self):
        self.assertEqual(jaccard_edges({}, .8), ([], None))
        self.assertEqual(jaccard_edges({"a": {"x"}}, .8), ([], None))
        for fraction in (-1, 2, float("nan")):
            with self.assertRaises(ValueError):
                jaccard_edges({}, fraction)
        with self.assertRaises(ValueError):
            jaccard_edges({"a": set()}, .8)

    def fixture(self):
        pathways = [dict(go_id="A", description="A", genes=["x", "y"]),
                    dict(go_id="B", description="B", genes=["z"])]
        metadata = dict(go_id="A", annotation_size=2, tested_size=2, ct_count=1, ct_ratio=.5)
        tables = {m: [dict(metadata, **{"hypergeometric_bh" if m == "hypergeometric" else "p_bh": .01 if m == "hypergeometric" else .8})] for m in METHODS}
        return pathways, tables

    def test_selection_and_historical_small_pathway(self):
        paths, tables = self.fixture()
        args = (paths, tables, {"x", "y", "z"}, {"x"}, "short")
        self.assertEqual([p["go_id"] for p in select_pathways(*args)[1]], ["A"])
        self.assertEqual(select_pathways(*args, method_rule="all")[1], [])
        audit, chosen = select_pathways(*args, requested={"B"})
        self.assertEqual([p["go_id"] for p in chosen], ["B"])
        self.assertEqual(audit[1]["sumstat_bh"], "")  # Unmeasured tests remain unavailable.
        self.assertEqual(audit[1]["selection_basis"], "provided_list")
        with self.assertRaises(ValueError):
            select_pathways(*args, requested={"unknown"})

    def test_rejects_mixed_metadata_and_bad_probabilities(self):
        paths, original = self.fixture()
        for field, value in (("ct_count", 2), ("ct_ratio", .1), ("p_bh", float("nan")), ("p_bh", 1.1)):
            tables = copy.deepcopy(original)
            tables["sumstat"][0][field] = value
            with self.assertRaises(ValueError):
                select_pathways(paths, tables, {"x", "y", "z"}, {"x"}, "short")
        original["gsea"] = []
        with self.assertRaises(ValueError):
            select_pathways(paths, original, {"x", "y", "z"}, {"x"}, "short")


if __name__ == "__main__":
    unittest.main()
