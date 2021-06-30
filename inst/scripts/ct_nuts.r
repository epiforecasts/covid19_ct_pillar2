options(echo = TRUE)

library("here")

source(here::here("R", "ct_model.r"))

suppressWarnings(dir.create(here::here("output")))

if (interactive()) {
  model <- "symptoms"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  model <- args[1]
}

for (target_gene in 1:2) {
  fit <- fit_model("symptoms", target_gene, threads = threading(36))
  saveRDS(fit, here::here("output", paste0("ct_", model, target_gene, ".rds")))
}
