options(echo = TRUE)

library("here")

suppressWarnings(dir.create(here::here("output")))

if (!interactive()) {
  source(here::here("R", "ct_model.r"))
  args <- commandArgs(trailingOnly = TRUE)
  model <- args[1]
}

for (target_gene in 1:2) {
  fit <- fit_model(model, target_gene, sample = 0.1, fit_to = "q", threads = threading(36))
  saveRDS(fit, here::here("output", paste0("ct_", model, target_gene, "_ct_subsampled.rds")))
}

for (target_gene in 1:2) {
  fit <- fit_model(model, target_gene, fit_to = "q", threads = threading(36))
  saveRDS(fit, here::here("output", paste0("ct_", model, target_gene, "_ct.rds")))
}
