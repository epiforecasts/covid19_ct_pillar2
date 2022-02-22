library(forcats)
library(dplyr)
library(socialmixr)
library(janitor)

load_data <- function(variant_file, max_days_since_onset = 6) {

  epidemic_phase <-
    c(Alpha = "Dec 2020",
      Delta = "Jun 2021",
      Omicron = "Dec 2021")
  declining_phase <-
    c(Alpha = "Jan 2021",
      Delta = "Dec 2021",
      Omicron = "Jan 2022")

  variants <- readRDS(variant_file) %>%
    mutate(
      variant = recode_factor(variant, Wildtype = "_Wildtype"),
      variant = fct_relevel(variant, "_Wildtype"),
      dose = factor(dose),
      days_since_onset = as.integer(date_specimen - onsetdate),
      dose_variant = paste(dose, variant, sep = "_"),
      reinfection = if_else(
        episode > 1, "reinfection", "no known previous infection"
      ),
      variant = sub("BA.2", "BA2", variant),
      ivdr = factor(interaction(variant, dose, reinfection)),
      sex = factor(sex),
      source_lab = factor(source_lab),
      month_year = format(onsetdate, "%b %Y"),
      epidemic = epidemic_phase[variant],
      declining = declining_phase[variant],
      phase = if_else(
        !is.na(epidemic) & epidemic == month_year, "increasing",
        if_else(!is.na(declining) & declining == month_year, "declining",
                NA_character_)),
      age_group = socialmixr::reduce_agegroups(
                                age, limits = seq(0, 80, by = 20)),
      age_group = socialmixr::limits_to_agegroups(age_group)
    ) %>%
    filter(days_since_onset >= 0 & days_since_onset <= max_days_since_onset)

  onset_week_start <- wday(max(variants$onsetdate) - 6, week_start = 1)
  variants <- variants %>%
    mutate(
      week_onset = floor_date(onsetdate, "week", onset_week_start)
    )

  variants_freq <- tabyl(variants, dose_variant) %>%
    filter(n > 1000) %>%
    select(dose_variant)

  variants <- variants %>%
    inner_join(variants_freq, by = "dose_variant") %>%
    select(source_lab, week_onset, days_since_onset, sex, age, age_group,
           reinfection, dose, variant, ivdr, phase, ct = p2ch2cq)
}
