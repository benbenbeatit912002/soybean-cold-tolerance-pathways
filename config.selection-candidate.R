config <- list(
  label="Unconfirmed candidate reconstruction of historical Excel selection",
  historical_procedure_confirmed=FALSE,
  count_method="hypergeometric", probability_field="p_value", alpha=.05, method_rule="any",
  phases=list(short=list(min_ct=6,min_ratio=.05,biological_rule="or"),
              middle=list(min_ct=40,min_ratio=0,biological_rule="count")))
