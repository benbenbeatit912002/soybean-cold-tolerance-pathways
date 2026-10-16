# A tiny teaching example. These are invented numbers, not research results.
# Run: Rscript examples/walkthrough.R
# Read examples/walkthrough.py beside this file to compare the same steps.
source("R/core.R")

# 1. Ten measured genes. A larger external score means stronger evidence here.
scores <- setNames(c(5, 0, 2, 7, 1, 4, 8, 3, 6, 2), LETTERS[1:10])

# 2. The CT list is separate from external score values.
ct_genes <- c("A", "D", "G")
pathways <- data.frame(go_id = "toy:1", description = "Example pathway",
                       genes = I(list(c("A", "D", "G"))))

# 3. Compute all overlap counts in the same measured universe.
sets <- prepare_membership(pathways, names(scores), ct_genes, 2, 8)
cat("Measured genes:", attr(sets, "universe_size"), "\n")
cat("Observed SUMSTAT:", sum(scores[sets$genes[[1]]]), "\n")  # 5 + 7 + 8 = 20
cat("Hypergeometric p-value:", hypergeometric(sets), "\n")

# 4. Compare 20 with sums from 999 random three-gene sets.
print(permutation_test(scores, sets$genes, "sumstat", B = 999, seed = 17))
cat("The Python estimate can differ slightly: each language generates different random draws.\n")
