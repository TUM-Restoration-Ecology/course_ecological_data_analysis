Analysis of Ecological Data <br> 3.2 Select and prepare traits
================
<b>Markus Bauer</b> <br>
<b>2024-07-15</b>

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
library(GIFT)
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
  here("data", "processed", "data_processed_traits.csv"),
  col_names = TRUE, col_types = cols(.default = "?")
  )

species <- data_species %>%
  semi_join(data_traits, by = "name") # only herbal species
traits <- data_traits %>%
  semi_join(species, by = "name") # only species which occur in the plots

bodysize <- c(10, 20, 30, 40, 50, NA, 70)
carnivory <- c(1, 1, 0, 1, 0,1, 0)
red <- c(1, 0, 0.5, 0, 0.2, 0, 1)
yellow <- c(0, 1, 0, 0, 0.3, 1, 0)
blue <- c(0, 0, 0.5,1, 0.5, 0, 0)
colors.fuzzy <- cbind(red, yellow, blue)
names <- paste("sp", 1:7, sep="")
tall <- as_tibble(cbind(names, bodysize, carnivory, colors.fuzzy))
tall
```

    ## # A tibble: 7 × 6
    ##   names bodysize carnivory red   yellow blue 
    ##   <chr> <chr>    <chr>     <chr> <chr>  <chr>
    ## 1 sp1   10       1         1     0      0    
    ## 2 sp2   20       1         0     1      0    
    ## 3 sp3   30       0         0.5   0      0.5  
    ## 4 sp4   40       1         0     0      1    
    ## 5 sp5   50       0         0.2   0.3    0.5  
    ## 6 sp6   <NA>     1         0     1      0    
    ## 7 sp7   70       0         1     0      0

``` r
rm(list = setdiff(ls(), c("species", "traits", "tall")))
```

# 2 Different trait types

We have categorical (binary, nominal, ordinal), continuous, proportional
or circular traits. Reload the traits table and define the column type
with - “col_logical()” = binary: yes/no data or presence-absence data -
“col_integer()” = whole numbers: counting data or rounded data -
“col_double()” = continous numbers: measured values - “col_character()”
= text - “col_factor()” or specify further “col_factor(levels,
ordered)”: ordinal traits are hierarchically and nominal traits not -
“col_date()” or specify further “col_date(”%Y-%m-%d”)“: phenological
traits are such circular traits

``` r
data_traits <- read_csv(
  here("data", "processed", "data_processed_traits.csv"),
  col_names = TRUE, col_types = cols(
    .default = "?",
    name = col_character()#,
    # go on here
    )
  )
```

Fuzzy coding: Different individaulas of the same species belong to
different levels of a factor

# 3 Missing traits

Identify missing traits.

``` r
traits %>%
  select(accepted_name, starts_with("trait_value")) %>%
  filter(is.na(trait_value_1.6.2) | is.na(trait_value_1.6.3))
```

    ## # A tibble: 2 × 3
    ##   accepted_name       trait_value_1.6.2 trait_value_1.6.3
    ##   <chr>                           <dbl>             <dbl>
    ## 1 Aconitum lycoctonum                 2                NA
    ## 2 Molinia arundinacea                 2                NA

Possibility 1: Remove the row.

Possibility 2: Use the sample average

Possibility 3: Measure yourself in the field

Possibility 4: Use trait values from phylogenetically close species.
(Recommended)

Download the trait values and search for closely related species:

``` r
data_gift <- GIFT_traits(
  trait_IDs = c("1.2.2", "1.6.3","3.7.1", "3.21.1", "4.5.1"),
  agreement = 0.66, bias_ref = FALSE, bias_deriv = FALSE
  )

(data <- data_gift %>%
  filter(str_detect(work_species, "Aconitum")))
```
