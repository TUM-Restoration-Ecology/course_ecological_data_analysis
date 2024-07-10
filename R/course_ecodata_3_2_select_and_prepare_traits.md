Analysis of Ecological Data <br> 3.2 Select and prepare traits
================
<b>Markus Bauer</b> <br>
<b>2024-07-10</b>

- [1 Preparation](#1-preparation)
  - [1.1 The example data set](#11-the-example-data-set)
  - [1.2 Load data](#12-load-data)
- [2 Different trait types](#2-different-trait-types)
- [3 Missing traits](#3-missing-traits)

**Markus Bauer**

Technichal University of Munich, TUM School of Life Sciences, Chair of
Restoration Ecology, Emil-Ramann-Straße 6, 85354 Freising, Germany

<markus1.bauer@tum.de>

<img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" alt="ORCID logo" width="20" height="20"/><https://orcid.org/0000-0001-5372-4174>

[![Google
Scholar](https://img.shields.io/badge/Google%20Scholar-4285F4?style=for-the-badge&logo=google-scholar&logoColor=white)](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/markus1bauer)

# 1 Preparation

Load the libraries. You can always install missing libraries by
`install.packages("tidyverse")` for the tidyverse packages.

``` r
library(here)
library(tidyverse)
rm(list = ls())
```

## 1.1 The example data set

A study in the Ammer valley before and after the weir in the [floodplain
Schnalz](https://www.openstreetmap.org/#map=16/47.7737/10.9615).

Normally there are alluvial grey alder forests (Alnetum incanae,
alliance [Alnion
incanae](https://floraveg.eu/vegetation/overview/Alnion%20incanae),
EUNIS habitat [T12](https://floraveg.eu/habitat/overview/T12))

Due to the dam along the Ammer and the weir in the Ammer the flood
regime changed and the vegetation could have developed to a maple-ash
forest (‘Edellaubholzwald’) (Adoxo-Aceretum, alliance
[Tilio-Acerion](https://floraveg.eu/vegetation/overview/Tilio-Acerion),
EUNIS habitat type [T1F](https://floraveg.eu/habitat/overview/T1F%3E))

The dataset is also available on GitHub (Bauer et al. 2018)

[![DOI:10.14471/2018.38.006](http://img.shields.io/badge/DOI-10.14471/2018.38.006-informational.svg)](https://doi.org/10.14471/2018.38.006)

<https://github.com/markus1bauer/2018_alluvial_forest_river_ammer/tree/main>

## 1.2 Load data

``` r
data_species <- read_csv(
  here("data", "raw", "example_course_alluvial_forest_mb", "data_raw_species.csv"),
  col_names = TRUE, col_types = cols(
    .default = "?",
    name = "f",
    layer = "f"
    )
  ) %>%
  select(-contains("Extra")) %>% # exclude some plots
  filter(layer == "h") %>% # filter only species of the herbal layer
  group_by(name) %>%
  mutate(
    sum = sum(c_across(IN1:AC6)),
    presence = if_else(sum > 0, 1, 0)
    ) %>%
  filter(presence == 1) %>% # exclude species which are not present in the remaining plots
  ungroup() %>%
  select(-sum, -presence, -layer)

data_traits <- read_csv(
  here("data", "raw", "example_course_alluvial_forest_mb", "data_raw_traits.csv"),
  col_names = TRUE, col_types = cols(.default = "?")
  ) %>%
  select(name, lifeform) %>%
  filter(!str_detect(lifeform, "wood"))

species <- data_species %>%
  semi_join(data_traits, by = "name") # only herbal species
traits <- data_traits %>%
  semi_join(species, by = "name") # only species which occur in the plots

rm(list = setdiff(ls(), c("species", "traits")))
```

# 2 Different trait types

# 3 Missing traits
