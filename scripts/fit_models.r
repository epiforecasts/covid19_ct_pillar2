library(forcats)
library(dplyr)
library(janitor)
library(socialmixr)
library(lubridate)
library(mgcv)
library(tidyr)

source(here::here("R", "data.r"))
source(here::here("R", "utils.r"))
source(here::here("R", "filenames.r"))

variants <- load_data(variant_file)

suppressWarnings(dir.create(here::here("output")))

## GAM fit
fit <- bam(ct ~
             0 +
             ivdr +
             s(age) +
             s(sex, bs = 're') +
             s(source_lab, bs = 're') +
             s(days_since_onset, by = ivdr, k = 7),
           data = variants, family = Gamma(link = "identity"))
saveRDS(fit, fit_file)

## means estimates
means <- estimate_means(variants, c(
  "days_since_onset", "reinfection", "dose", "variant")
  ) %>%
  mutate(days_since_onset = as.integer(sub("^var", "", days_since_onset)))
age_means <- estimate_means(variants, c(
  "days_since_onset", "reinfection", "dose", "variant", "age_group")
  ) %>%
  mutate(days_since_onset = as.integer(sub("^var", "", days_since_onset)))
phase_means <-
  estimate_means(variants, c(
    "days_since_onset", "reinfection", "dose", "variant", "phase")
    ) %>%
  mutate(days_since_onset = as.integer(sub("^var", "", days_since_onset)))
age_phase_means <-
  estimate_means(variants, c(
    "days_since_onset", "reinfection", "dose", "variant", "phase", "age_group")
    ) %>%
  mutate(days_since_onset = as.integer(sub("^var", "", days_since_onset)))
time_means <-
  estimate_means(variants, c(
    "week_onset", "variant")
    )
saveRDS(list(means = means,
             age_means = age_means,
             phase_means = phase_means,
             age_phase_means = age_phase_means,
             time_means = time_means),
        means_file)
