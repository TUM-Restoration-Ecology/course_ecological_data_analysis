Analysis of Ecological Data <br> 3.1 How to get traits?
================
<b>Markus Bauer</b> <br>
<b>2024-07-10</b>

- [1 Preparation](#1-preparation)
  - [1.2 The example data set](#12-the-example-data-set)
  - [1.2 Load data](#12-load-data)
- [2 Metadata of plant traits of
  GIFT](#2-metadata-of-plant-traits-of-gift)
- [3 Use TNRS for name resolving](#3-use-tnrs-for-name-resolving)
- [4 Download trait values](#4-download-trait-values)
  - [4.1 Species level](#41-species-level)
  - [4.2 Taxonomic level](#42-taxonomic-level)
- [5 Merge your dataset with GIFT
  traits](#5-merge-your-dataset-with-gift-traits)
  - [5.1 Merge dataset](#51-merge-dataset)
  - [5.1 Check completeness](#51-check-completeness)
- [6 Save your traits table](#6-save-your-traits-table)

<br/> <br/> <b>Markus Bauer</b>

Technichal University of Munich, TUM School of Life Sciences, Chair of
Restoration Ecology, Emil-Ramann-Straße 6, 85354 Freising, Germany

<markus1.bauer@tum.de>

<a href="https://orcid.org/0000-0001-5372-4174"><img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" alt="ORCID logo" width="20" height="20"/>
<https://orcid.org/0000-0001-5372-4174>

[![Google
Scholar](https://img.shields.io/badge/Google%20Scholar-4285F4?style=for-the-badge&logo=google-scholar&logoColor=white)](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/markus1bauer)

# 1 Preparation

Load the libraries. You can always install missing libraries by
`install.packages("GIFT")` for the GIFT package.

``` r
library(here)
library(tidyverse)
library(GIFT)
library(TNRS)
rm(list = ls())
```

## 1.2 The example data set

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

The dataset is also available on GitHub Bauer et al. (2018) Tuexenia
<br>
[![DOI:10.14471/2018.38.006](http://img.shields.io/badge/DOI-10.14471/2018.38.006-informational.svg)](https://doi.org/10.14471/2018.38.006):
<br>
<https://github.com/markus1bauer/2018_alluvial_forest_river_ammer/tree/main>

## 1.2 Load data

Load the example dataset of alluvial forests along river Ammer.

``` r
species <- read_csv(
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

traits <- read_csv(
  here("data", "raw", "example_course_alluvial_forest_mb", "data_raw_traits.csv"),
  col_names = TRUE, col_types = cols(.default = "?")
  ) %>%
  select(name, lifeform) %>%
  filter(!str_detect(lifeform, "wood"))

species <- species %>%
  semi_join(traits, by = "name") # only herbal species
traits <- traits %>%
  semi_join(species, by = "name") # only species which occur in the plots
```

# 2 Metadata of plant traits of GIFT

There are many functional traits available in GIFT. Each of these traits
has an identification number called `trait_ID`. Since the two functions
for retrieving trait values, `GIFT_traits()` and `GIFT_traits_raw()`,
rely on these IDs, the first step is to call the function
`GIFT_traits_meta()` to know what the ID of the desired trait is. For
example, let’s say we want to retrieve the maximum vegetative heights of
plant species.

``` r
trait_meta <- GIFT::GIFT_traits_meta()
trait_meta %>%
  filter(str_detect(Trait2, "height"))
```

    ##   Lvl1     Category Lvl2       Trait1   Lvl3            Trait2 Units    type
    ## 1    1   Morphology  1.6 Plant height  1.6.1  Plant_height_min     m numeric
    ## 2    1   Morphology  1.6 Plant height  1.6.2  Plant_height_max     m numeric
    ## 3    1   Morphology  1.6 Plant height  1.6.3 Plant_height_mean     m numeric
    ## 4    3 Reproduction 3.12  Seed height 3.12.1   Seed_height_min    mm numeric
    ## 5    3 Reproduction 3.12  Seed height 3.12.2   Seed_height_max    mm numeric
    ## 6    3 Reproduction 3.15 Fruit height 3.15.1  Fruit_height_min    cm numeric
    ## 7    3 Reproduction 3.15 Fruit height 3.15.2  Fruit_height_max    cm numeric
    ## 8    3 Reproduction 3.15 Fruit height 3.15.3 Fruit_height_mean    cm numeric
    ##   comment count
    ## 1    <NA> 43353
    ## 2    <NA> 70367
    ## 3    <NA> 23219
    ## 4    <NA>     7
    ## 5    <NA>   148
    ## 6    <NA>    54
    ## 7    <NA>    54
    ## 8    <NA>    54

We can see that the ID of this trait is 1.6.2. Now that we have the ID,
we can use `GIFT_traits()` to retrieve the growth form values for
different plant species.

Search the IDs of specific leaf area (SLA) and seed mass.  
Note the categories, IDs and unit for the material-and-methods section.

# 3 Use TNRS for name resolving

``` r
metadata <- TNRS::TNRS_metadata()

names <- species %>%
  rowid_to_column("id") %>%
  select(id, name) %>%
  TNRS::TNRS()
names2 <- names %>%
  select(
    Name_submitted, Taxonomic_status, Accepted_name, Accepted_name_url,
    Accepted_family
    ) %>%
  rename_with(tolower)
```

Note the version of the database from `metadata$version` and the sources
of the names with `metadata$citations` for the material-and-methods
section

Note the range of `Overall_score` for your submitted species names.

# 4 Download trait values

## 4.1 Species level

## 4.2 Taxonomic level

# 5 Merge your dataset with GIFT traits

## 5.1 Merge dataset

## 5.1 Check completeness

# 6 Save your traits table
