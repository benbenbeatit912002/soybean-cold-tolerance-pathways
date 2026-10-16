"""A tiny teaching example. These are invented numbers, not research results.

Run: python examples/walkthrough.py
Then read the corresponding examples/walkthrough.R side by side.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"python"))
from ctpathway import prepare_membership, hypergeometric, permutation_test

# 1. Ten measured genes. A larger external score means stronger evidence here.
scores = dict(zip("ABCDEFGHIJ", [5, 0, 2, 7, 1, 4, 8, 3, 6, 2]))

# 2. A CT list is a membership list, independent of the external score values.
ct_genes = ["A", "D", "G"]
pathways = [{"go_id": "toy:1", "description": "Example pathway", "genes": ["A", "D", "G"]}]

# 3. Use the same ten-gene universe for pathway size and CTgene counts.
sets, universe_size, ct_count = prepare_membership(pathways, scores, ct_genes, 2, 8)
print("Measured genes:", universe_size)
print("Observed SUMSTAT:", sum(scores[g] for g in sets[0]["genes"]))  # 5 + 7 + 8 = 20
print("Hypergeometric p-value:", hypergeometric(sets, universe_size, ct_count)[0])

# 4. Compare 20 with sums from 999 random three-gene sets.
result = permutation_test(scores, [sets[0]["genes"]], "sumstat", permutations=999, seed=17)
print("SUMSTAT random-set result:", result[0])
print("The R estimate can differ slightly: each language generates different random draws.")
