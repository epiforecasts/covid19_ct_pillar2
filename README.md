# Cycle threshold values in symptomatic Covid-19 cases in England

In this work we analyse Ct values as a proxy for viral loads in symptomatic Pillar 2 cases in England. We describe trends as a function of the number of days since symptom onset, age, variant, number of vaccine doses, reinfection status, processing laboratory, and interactions between these variables. We further try to recover time trends using a model that incorporates all these factors, and investigate the role of epidemic phase bias.

This is work in progress. Suggestions, comments or corrections are be very welcome.

The report is available in [html](https://epiforecasts.io/covid19_ct_pillar2/report) or [pdf](https://epiforecasts.io/covid19_ct_pillar2/report.pdf) format.

The analysis is implemented in `R` with heavy use of the `mgcv` package.

# Citation

S. Funk, S. Abbott. _Cycle threshold values in symptomatic Covid-19 cases in England_ (2022). https://github.com/epiforecasts/covid19_ct_pillar2

# Documentation

## File structure

Folder | Purpose
---|---
[`R`](R/) | R functions for preprocessing and analysing data.
[`scripts`](scripts/) | Scripts used to create results.

## Obtaining data

The Pillar 2 surveillance data used in this repository is currently not publicly available. Access has been provided by the UK Health Security Agency (UKHSA) through the Scientific Pandemic Influenza Group on Modelling, Operational sub-group (SPI-M-O) for the Scientific Advisory Group for Emergencies (SAGE), where this work was first presented.

## Dependencies

All dependencies can be installed using the following,

```{r}
remotes::install_dev_deps()
```

## Reproducing our results

Once data access has been obtained run `scripts/fit_model.r`. Note that fitting the model may take up to 10 hours on a fast multi-core machine. Alternatively, render `docs/report.Rmd` and this script will be called if results are not present.
