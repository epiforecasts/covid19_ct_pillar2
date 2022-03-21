# Viral loads in symptomatic Covid-19 cases in England

In this work we analyse Ct values as a proxy for viral loads in symptomatic Pillar 2 cases in England. We describe trends as a function of the number of days since symptom onset, age, variant, number of vaccine doses, reinfection status, and processing laboratory. We further try to recover time trends using a model that incorporates all these factors, and investigate the role of epidemic phase bias.

The analysis is implemented in `R` with heavy use of the `mgcv` package.

# Citation

S. Funk, S. Abbott. _Viral loads in symptomatic Covid-19 cases in England _ (2022). https://github.com/epiforecasts/covid19_ct_pillar2

# Documentation

## File structure

Folder | Purpose
---|---
[`R`](R/) | R functions for preprocessing and analysing data.
[`scripts`](scripts/) | Scripts used to create results.

## Obtaining data

The data used in this repository is currently not publicly available.

## Dependencies

All dependencies can be installed using the following, 

```{r}
remotes::install_dev_deps()
```
