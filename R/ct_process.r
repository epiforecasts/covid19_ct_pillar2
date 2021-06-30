library("dplyr")
library("tidyr")
library("tidybayes")
library("ggplot2")
library("brms")
library("loo")

suppressWarnings(dir.create(here::here("figure")))

variants <- readRDS(here::here("data", "processed", "variants_vl.rds")) %>%
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

output_files <- list.files(here::here("output"), pattern = "^ct_") %>%
  sub(pattern = "ct_", replacement = "") %>%
  sub(pattern = "[0-9]\\.rds", replacement = "") %>%
  unique()

loos <- list()
draws <- list()
ages <- list()
effects <- list()

for (file in output_files) {
  loos[[file]] <- list()
  draws[[file]] <- list()
  ages[[file]] <- list()
  file_names <- here::here("output", paste0("ct_", file, 1:2, ".rds"))
  names(file_names) <- c("ORF1ab", "N")
  for (gene in names(file_names)) {
    fit_file <- file_names[gene]
    if (file.exists(fit_file)) {
      fit <- readRDS(fit_file)
      draws[[file]][[gene]] <- fit %>%
        recover_types(variants) %>%
        spread_draws(`b_.*`, regex = TRUE) %>%
        median_qi(.width = c(.90, .50))
      if ("days_since_onset" %in% names(fit$data)) {
        fit$data$days_since_onset <- as.integer(fit$data$days_since_onset )
      }
      ages[[file]][[gene]] <- conditional_smooths(fit) %>%
        .$mu %>%
        tibble()
    }
  }
  draws[[file]] <- bind_rows(draws[[file]], .id = "target_gene")
  ages[[file]]<- bind_rows(ages[[file]], .id = "target_gene")
  effects[[file]] <- draws[[file]] %>%
    pivot_longer(starts_with("b_")) %>%
    filter(!grepl("_Intercept", name)) %>%
    mutate(name = sub("\\.(lower|upper)$", "|\\1", name),
           name = if_else(grepl("\\|", name), name,
                          paste(name, .point, sep = "|"))) %>%
    separate(name, c("variable", "name"), sep = "\\|") %>%
    pivot_wider() %>%
    mutate(variable = sub("^b_", "", variable),
           variable = gsub("_", " ", variable),
           variable = if_else(grepl("^dose ", variable),
                              sub("dose vaccinedose", "Dose ", variable),
                              variable),
           variable = if_else(grepl("Kent", variable), "Variant Alpha",
                              variable),
           variable = if_else(grepl("India", variable), "Variant Delta",
                              variable),
           variable = if_else(variable == "sexMale", "Sex: male", variable),
           variable = if_else(variable == "Rt", "R", variable),
           variable = if_else(variable == "days since onset", "Days since onset", variable),
           variable = if_else(variable == "cases", "Cases", variable)) %>%
    dplyr::mutate(type = dplyr::if_else(grepl("^Dose", variable),
                                        "Vaccine", "Epidemiology"),
                  type = dplyr::if_else(grepl("variant$", variable),
                                        "Variant", variable),
                  type = dplyr::if_else(grepl("variant$", variable),
                                        "Variant", variable),
                  type = dplyr::if_else(variable == "Male",
                                        "Sex", variable),
                  ) %>%
    dplyr::arrange(type, variable) %>%
    dplyr::mutate(variable = factor(variable, levels = unique(variable)),
                  variable = forcats::fct_rev(variable))
  p <- ggplot2::ggplot(effects[[file]], aes(x = variable, y = median,
                               ymin = lower, ymax = upper, color = target_gene)) +
    tidybayes::geom_pointinterval(position = position_dodge(width = 0.5)) +
    ggplot2::scale_colour_brewer("Target gene", palette = "Dark2") +
    ggplot2::xlab("") +
    ggplot2::ylab("Estimated additive effect on log10 viral load") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  ggsave(here::here("figure", paste0("ct_effects_", file, ".pdf")), p,
         width = 5, height = 7)

  p <- ggplot2::ggplot(ages[[file]], aes(x = age, y = estimate__,
                                    ymin = lower__, ymax = upper__,
                                    colour = target_gene, fill = target_gene)) +
    ggplot2::scale_colour_brewer("Target gene", palette = "Dark2") +
    ggplot2::scale_fill_brewer("Target gene", palette = "Dark2") +
    ggplot2::xlab("Age") +
    ggplot2::ylab("Estimated additive effect on log10 viral load") +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(alpha = 0.35, colour = NA) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  ggsave(here::here("figure", paste0("ct_age_", file, ".pdf")), p,
         width = 5, height = 5)

}

