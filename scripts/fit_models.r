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

ct <- load_data(ct_file)

suppressWarnings(dir.create(here::here("output")))

## GAM fit
system.time(fit <- bam(ct ~
             0 +
             variant +
             dose +
	     reinfection +
             sex +
             source_lab +
	     s(age, k = 10) +
             s(days_since_onset, k = 7),
             data = ct, family = Gamma(link = "identity")))
saveRDS(fit, fit_file)

## means estimates
means <- estimate_means(ct, c(
  "days_since_onset", "reinfection", "dose", "variant")
  ) %>%
  mutate(days_since_onset = as.integer(sub("^var", "", days_since_onset)))
age_means <- estimate_means(ct, c(
  "days_since_onset", "reinfection", "dose", "variant", "age_group")
  ) %>%
  mutate(days_since_onset = as.integer(sub("^var", "", days_since_onset)))
phase_means <-
  estimate_means(ct, c(
    "days_since_onset", "reinfection", "dose", "variant", "phase")
    ) %>%
  mutate(days_since_onset = as.integer(sub("^var", "", days_since_onset)))
age_phase_means <-
  estimate_means(ct, c(
    "days_since_onset", "reinfection", "dose", "variant", "phase", "age_group")
    ) %>%
  mutate(days_since_onset = as.integer(sub("^var", "", days_since_onset)))
time_means <-
  estimate_means(ct, c(
    "week_onset", "variant")
    )
saveRDS(list(means = means,
             age_means = age_means,
             phase_means = phase_means,
             age_phase_means = age_phase_means,
             time_means = time_means),
        means_file)
