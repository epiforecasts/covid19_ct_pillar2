estimate_means <- function(x, by) {

  x[["var"]] <- do.call(interaction, lapply(by, function(y) x[[y]]))

  means_fit <- mgcv::bam(ct ~ 0 + var, data = filter(x, !is.na(var)))
  ## get CI
  means <- coef(means_fit)
  V <- vcov(means_fit, unconditional = TRUE)
  se <- sqrt(diag(V))

  tb <- tibble(
    var = names(means),
    estimate = means,
    se = se,
    lower = means - 2 * se,
    upper = means + 2 * se
  )

  tb %>%
    separate(var, by, sep = "\\.")
}
