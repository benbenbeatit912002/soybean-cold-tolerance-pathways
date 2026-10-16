"""Independent small examples and edge cases for the Python implementation."""
import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
import numpy as np
from ctpathway import (prepare_membership, hypergeometric, bh_adjust, transform_score,
                      monte_carlo_p, gsea_es, maxmean_stat, fosco, fit_fosco,
                      resolve_duplicate_scores, permutation_test)


class StatisticsTests(unittest.TestCase):
    def test_hypergeometric_by_enumeration(self):
        genes = list("ABCDEFGHIJ")
        pathways = [{"go_id": "GO:1", "description": "example", "genes": list("ABDX")}]
        sets, n, k = prepare_membership(pathways, genes, list("ABC"), 2, 4)
        draws = list(itertools.combinations(genes, 3))
        exact = sum(len(set(draw) & set("ABD")) >= 2 for draw in draws)/len(draws)
        self.assertAlmostEqual(hypergeometric(sets, n, k)[0], exact)
        self.assertEqual(sets[0]["annotation_size"], 4)
        self.assertEqual(sets[0]["tested_size"], 3)

    def test_bh_known_example(self):
        np.testing.assert_allclose(bh_adjust([.01, .04, .03, .002]), [.02, .04, .04, .008])

    def test_score_domains(self):
        np.testing.assert_allclose(transform_score([.1, 1], "pvalue"), [1, 0])
        np.testing.assert_allclose(transform_score([2, .5], "fold_change"), [np.log(2)]*2)
        for values, kind in [([0, .5], "pvalue"), ([1, np.inf], "score"), ([-1], "score")]:
            with self.assertRaises(ValueError):
                transform_score(values, kind)

    def test_fosco_stable_and_domains(self):
        np.testing.assert_allclose(fosco([.1, .5], [4, 9], .5), 1-(1-np.array([.1, .5]))**[2, 3])
        self.assertGreater(fosco([1e-30], [2], .5)[0], 0)
        self.assertEqual(fosco([1], [2], .5)[0], 1)
        with self.assertRaises(ValueError):
            fosco([.1], [0], .5)
        with self.assertRaises(ValueError):
            fit_fosco([.1, .2], [2, 2])

    def test_gsea_against_full_running_sum(self):
        scores = np.array([5, 4, 3, 2, 1, 0])
        for hits in itertools.combinations(range(len(scores)), 2):
            increments = np.full(len(scores), -1/(len(scores)-len(hits)))
            increments[list(hits)] = scores[list(hits)]/sum(scores[list(hits)])
            self.assertAlmostEqual(gsea_es(scores, hits), max(0, np.cumsum(increments).max()))
        self.assertEqual(gsea_es(np.zeros(6), [0, 1]), 1)
        with self.assertRaises(ValueError):
            gsea_es(scores, range(6))

    def test_ties_and_duplicates(self):
        self.assertEqual(monte_carlo_p(5, [1, 2, 3, 4]), 1/5)
        self.assertEqual(monte_carlo_p(1, [1]*4), 1)
        self.assertEqual(maxmean_stat([-2, 1, 0, 3]), 1)
        with self.assertRaises(ValueError):
            resolve_duplicate_scores(["a", "a"], [1, 3])
        self.assertEqual(resolve_duplicate_scores(["a", "a"], [1, 3], "mean"), {"a": 2})

    def test_seeded_permutations(self):
        scores = dict(zip("ABCDEFGHIJ", range(1, 11)))
        members = [list("ABD"), list("EF")]
        for method in ("sumstat", "gsea", "maxmean"):
            first = permutation_test(scores, members, method, 39, 17)
            second = permutation_test(scores, members, method, 39, 17)
            self.assertEqual(first, second)
            for row in first:
                self.assertTrue(1/40 <= row["p_value"] <= row["p_bh"] <= 1)
        with self.assertRaises(ValueError):
            permutation_test(scores, members, "sumstat", 0)


if __name__ == "__main__":
    unittest.main()
