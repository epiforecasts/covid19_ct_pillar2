options(echo = TRUE)

library("brms")
library("dplyr")
library("lubridate")
library("cmdstanr")
library("janitor")

models <-
  list(baseline = bf(log10vl ~ variant_alt_name + sex + s(age),
                     family = student(),
                     sigma ~ (1 | age_group)),
       symptoms = bf(log10vl ~
                       sex + age_group + dose_vaccine_variant + (days_since_onset):(sex + age_group + dose_vaccine_variant),
                     family = student(),
                     sigma ~ (1 | age_group)),
       symptoms_smooth = bf(log10vl ~
                              s(age) + sex + dose_vaccine_variant + days_since_onset:(sex + age_group + dose_vaccine_variant),
                     family = student(),
                     sigma ~ s(age)),
       symptoms_random = bf(log10vl ~
                        sex + days_since_onset + (1 + days_since_onset | age_group) + (1 + days_since_onset | dose_vaccine_variant) + sex:days_since_onset,
                        family = student(),
                        sigma ~ (1 | age_group)),
       symptoms_random_age = bf(log10vl ~ sex + days_since_onset + (1 + days_since_onset | variant_alt_name:dose:age_group) + sex:days_since_onset,
                                family = student(),
                                sigma ~ (1 | age_group)))

variants <- readRDS(here::here("data", "variants_vl.rds")) %>%
  filter(sex != "Unknown") %>%
  mutate(variant_alt_name = recode_factor(variant_alt_name,
                                          Wildtype = "_Wildtype"),
         dose_vaccine =
           if_else(dose == "dose0", "dose0",
                   paste(dose, product_display_type, sep = "_")),
         dose_vaccine = factor(dose_vaccine),
         onsetdate = dmy(onsetdate),
         days_since_onset = as.integer(date_specimen - onsetdate)) %>%
  filter(days_since_onset >= 0 & days_since_onset <= 28) %>%
  mutate(dose_vaccine_variant =
           paste(dose_vaccine, variant_alt_name, sep = "_"))

variants_freq <- tabyl(variants, dose_vaccine_variant) %>%
  filter(n > 1000) %>%
  select(dose_vaccine_variant)

variants <- variants %>%
  inner_join(variants_freq, by = "dose_vaccine_variant")

fit_model <- function(model, target_gene, fit_to = "vl", sample = 1, ...) {
  if (sample == 1) {
	  data <- variants
  } else {
	  data <- sample_frac(variants, sample)
  }
  fit <- brm(models[[model]],
             data = data %>%
               rename(log10vl =
                        !!quo_name(paste0("p2ch", target_gene, "c", fit_to))),
             chains = 2,
             cores = 2,
             iter = 2000,
             backend = "cmdstanr",
             ...)
  return(fit)
}

