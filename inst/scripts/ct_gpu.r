options(echo = TRUE)

library("here")

source(here::here("R", "ct_model.r"))

suppressWarnings(dir.create(here::here("output")))

for (target_gene in 1:2) {
  data <- variants %>%
     rename(log10vl = !!quo_name(paste0("p2ch", target_gene, "cvl")))
  model_code <- make_stancode(models[["symptoms"]], data)
  model_data <- make_standata(models[["symptoms"]], data)
  class(model_data) <- "list"
  model_file <- write_stan_file(model_code)
  mod_cl <- cmdstan_model(model_file, cpp_options = list(stan_opencl = TRUE))
  fit_cl <- mod_cl$sample(data = model_data,
                          chains = 2, parallel_chains = 2,
                          opencl_ids = c(0, 0), refresh = 0)
  saveRDS(fit_cl, here::here("output", paste0("ct_", model, target_gene, "_gpu.rds")))
}

