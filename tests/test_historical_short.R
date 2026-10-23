# Base R checks for the independently implemented historical calculator.
source("R/historical_short.R")
rows <- data.frame(source_row = 2:5, gene_id = paste0("g", 1:4),
                   gene_length = c(1, 2000, 2001, 98001), snp_markers = c(4, 0, 0, 0))
imputed <- historical_snp_imputation(rows)
stopifnot(identical(imputed$values, c(4, 2, 0, 0)),
          imputed$audit$mean_including_zeros[1] == 2,
          abs(historical_fosco(.2, 4, .5)-.36) < 1e-15,
          historical_fosco(.2, 0, .5) == 0,
          historical_fosco(1.14e-14, 10, .01) != -expm1(10^.01*log1p(-1.14e-14)))
export <- historical_score_export(c(1, 0, .9947477, 1.176836406102666e-14))
stopifnot(identical(export$display, c("1.00E-16", "Inf", "0.002287056", "13.9292839")))
fit <- fit_historical_short(c(.2, .4, .1, .8, .6), c(1, 8, 4, 0, 3))
stopifnot(nrow(fit$grid) == 100L, fit$lambda %in% ((1:100)/100),
          is.infinite(fit$score[4]), !anyNA(fit$adjusted))
failed <- tryCatch({ fit_historical_short(c(.2, .3), c(1, -1)); FALSE }, error = function(e) TRUE)
stopifnot(failed)
cat("Historical short recalculation tests passed.\n")
