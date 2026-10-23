# Synthetic pipeline demonstration

These files contain 20 invented genes named TOY01 through TOY20 and three invented pathways. They are not soybean measurements and must not be used as research results.

Short scores increase from 1 to 20; middle scores decrease from 20 to 1. Short CT membership is the last four genes and middle membership the first four. Prior scores are placeholders equal to one and are not added to the analysis scores.

TOY:1 comprises the last six genes; its short SUMSTAT is 15+16+17+18+19+20 = 105. TOY:2 comprises the first six genes. TOY:3 contains every second gene starting with TOY01.

Run `python run_pipeline.py --mode demo` or `Rscript --vanilla run_pipeline.R --mode demo` from the repository root. Both use 199 draws per distinct size/method, two phases, and the same output format as real-data analyses.
