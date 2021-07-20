options(echo = TRUE)

library("here")

if (!interactive()) {
  source(here::here("R", "ct_model.r"))
  args <- commandArgs(trailingOnly = TRUE)
  model <- args[1]
}

suppressWarnings(dir.create(here::here("output")))

for (target_gene in 1:2) {
  fit <- fit_model(model, target_gene, fit_to = "q", algorithm = "meanfield")
  saveRDS(fit, here::here("output", paste0("ct_", model, target_gene, "_ct_vb.rds")))
}
