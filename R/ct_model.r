options(echo = TRUE)

library("brms")
library("dplyr")
library("lubridate")
library("cmdstanr")

models <-
  list(baseline = bf(log10vl ~ variant_alt_name + sex + s(age),
                     family = student(),
                     sigma ~ (1 | age_group)),
       symptoms = bf(log10vl ~
                       variant_alt_name + sex + s(age) + days_since_onset:variant_alt_name + dose_vaccine_variant,
                     family = student(),
                     sigma ~ (1 | age_group)))

variants <- readRDS(here::here("data", "variants_vl.rds")) %>%
  filter(is.na(product_display_type) | product_display_type != "MD",
         sex != "Unknown") %>%
  rename(Rt = median) %>%
  mutate(dose_vaccine =
           if_else(dose == "dose0", "dose0",
                   paste(dose, product_display_type, sep = "_")),
         dose_vaccine = factor(dose_vaccine),
         onsetdate = dmy(onsetdate),
         days_since_onset = date_specimen - onsetdate) %>%
  filter(days_since_onset >= 0 & days_since_onset <= 28, dose != "dose2",
         !(variant_alt_name == "Wildtype" & dose != "dose0")) %>%
  mutate(dose_vaccine_variant =
           if_else(dose_vaccine == "dose0", "_none",
                   paste(dose_vaccine, variant_alt_name, sep = "_")))

fit_model <- function(model, target_gene, ...) {
  fit <- brm(models[[model]],
             data = variants %>%
               rename(log10vl = !!quo_name(paste0("p2ch", target_gene, "cvl"))),
             chains = 2,
             cores = 2,
             iter = 2000,
             backend = "cmdstanr",
             ...)
  return(fit)
}

