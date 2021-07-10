options(echo = TRUE)

library("here")

suppressWarnings(dir.create(here::here("output")))

if (!interactive()) {
  source(here::here("R", "ct_model.r"))
  args <- commandArgs(trailingOnly = TRUE)
  model <- args[1]
}

for (target_gene in 1:2) {
  fit <- fit_model(model, target_gene, sample = 0.1, threads = threading(36))
  saveRDS(fit, here::here("output", paste0("ct_", model, "_subsampled", target_gene, ".rds")))
}

for (target_gene in 1:2) {
  fit <- fit_model(model, target_gene, threads = threading(36))
  saveRDS(fit, here::here("output", paste0("ct_", model, target_gene, ".rds")))
}
