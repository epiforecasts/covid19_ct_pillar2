options(echo = TRUE)

library("here")

source(here::here("R", "ct_model.r"))

suppressWarnings(dir.create(here::here("output")))

for (target_gene in 1:2) {
  fit <- fit_model("symptoms", target_gene, method = "meanfield")
  saveRDS(fit, here::here("output", paste0("ct_", model, target_gene, "_vb.rds")))
}
