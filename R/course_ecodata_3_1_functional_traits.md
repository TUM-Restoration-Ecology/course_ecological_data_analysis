Analysis of Ecological Data <br> 3.1 How to get traits?
================
<b>Markus Bauer</b> <br>
<b>2024-07-10</b>

- [1 Preparation](#1-preparation)
  - [1.1 The example data set](#11-the-example-data-set)
  - [1.2 Load data](#12-load-data)
- [2 Metadata of traits](#2-metadata-of-traits)
- [3 Get trait values](#3-get-trait-values)
  - [3.1 How to download traits](#31-how-to-download-traits)
  - [3.2 Get the references of the
    measurements](#32-get-the-references-of-the-measurements)
  - [3.3 Download traits of SLA, height, seed
    mass](#33-download-traits-of-sla-height-seed-mass)
- [4 Name resolving](#4-name-resolving)
- [5 Merge](#5-merge)
- [6 Check completeness](#6-check-completeness)
- [7 Save](#7-save)

**Markus Bauer**

Technichal University of Munich, TUM School of Life Sciences, Chair of
Restoration Ecology, Emil-Ramann-Straße 6, 85354 Freising, Germany

<markus1.bauer@tum.de>

<img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" alt="ORCID logo" width="20" height="20"/><https://orcid.org/0000-0001-5372-4174>

[![GoogleScholar](https://img.shields.io/badge/Google%20Scholar-4285F4?style=for-the-badge&logo=google-scholar&logoColor=white)](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/markus1bauer)

# 1 Preparation

Load the libraries. You can always install missing libraries by
`install.packages("GIFT")` for the GIFT package.

``` r
library(here)
library(tidyverse)
library(GIFT)
library(TNRS)
library(naniar)
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

# 2 Metadata of traits

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

- Search the IDs of *mean* specific leaf area (SLA), seed mass aand
  plant height.

- Note the categories, IDs and unit for the material-and-methods
  section.

# 3 Get trait values

[Original
vignette](https://biogeomacro.github.io/GIFT/articles/GIFT.html#trait-data)
(Denelle & Weigelt 2024)

Get the functional trait values from the database GIFT ([Weigelt et
al. 2023](https://doi.org/10.1111/jbi.13623)) with the package GIFT
([Denelle et al. 2023](https://doi.org/10.1111/2041-210X.14213)).

## 3.1 How to download traits

There are two functions to access trait values. First,
`GIFT_traits_raw()` returns all trait values for a given species and a
given trait. These trait values can then vary. Second, `GIFT_traits()`
returns an aggregated trait value at the species level. The aggregation
simply takes the mean for continuous traits or the most frequent entry
for categorical traits. However, for some specific cases, the
aggregation takes either the minimum or the maximum, as for the trait we
chose. Let’s retrieve the raw and aggregated values for the maximum
vegetative height of plants (trait_ID 1.6.2).

``` r
height <- GIFT_traits(
  trait_IDs = c("1.6.2"), agreement = 0.66, bias_ref = FALSE, bias_deriv = FALSE
  )

height_raw <- GIFT_traits_raw(trait_IDs = c("1.6.2"))

# Raw values
height_raw %>%
  select(work_species, trait_ID, trait_value) %>%
  filter(work_species == "Fagus sylvatica")

# Aggregated value
height %>%
  filter(work_species == "Fagus sylvatica")
```

There were three maximum heights for Fagus sylvatica, 30, 35, and 50
meters, which led to an aggregated value of 50 meters.

## 3.2 Get the references of the measurements

If you want to look up the references that led to the aggregated trait
value, you can run this chunk:

``` r
references <- GIFT_references(GIFT_version = "latest")

unique(unlist(strsplit(height$references_1.6.2, ",")))

references <- references[
  which(references$ref_ID %in% 
          unique(unlist(strsplit(height$references_1.6.2, ",")))), ]
references[1:2, ]
```

## 3.3 Download traits of SLA, height, seed mass

- Get the aggregated trait values of the mean SLA, plant height and seed
  mass.

Here is an example how to download several traits:

``` r
data_gift <- GIFT_traits(
  trait_IDs = c("1.6.2", "1.6.3"), agreement = 0.66, bias_ref = FALSE, bias_deriv = FALSE
  )

rm(list = setdiff(ls(), c("species", "traits", "data_gift")))
```

# 4 Name resolving

[Original
vignette](https://cran.r-project.org/web/packages/TNRS/vignettes/TNRS_vignette.html)
(Maitner 2024)

Use the package TNRS ([Boyle et
al.2013](https://doi.org/10.1186/1471-2105-14-16)) for name resolving.

``` r
metadata <- TNRS::TNRS_metadata()

data <- species %>%
  rowid_to_column("id") %>%
  select(id, name) %>%
  TNRS::TNRS(
    sources = c("wcvp", "wfo"),
    classification = "wfo",
    mode = "resolve"
  )
names <- data %>%
  select(
    Name_submitted, Taxonomic_status, Accepted_name, Accepted_name_url,
    Accepted_family
    ) %>%
  rename_with(tolower) %>%
  rename(name = name_submitted)
```

- Note the version of the database from metadata\$version and the
  sources of the names with metadata\$citations for the
  material-and-methods section

- Note the range of Overall_score for your submitted species names.

# 5 Merge

Merge the resolved species names and the functional plant traits into
the existing tables ‘species’ and ‘traits’.

``` r
data <- species %>%
  full_join(names %>% select(name, accepted_name), by = "name") %>%
  select(-name)
species <- data

data <- names %>%
  full_join(traits, by = "name") %>%
  left_join(data_gift %>% rename(accepted_name = work_species), by = "accepted_name")
traits <- data %>%
  select(-work_ID, -work_author, -starts_with("references"))

rm(list = setdiff(ls(), c("species", "traits")))
```

# 6 Check completeness

``` r
traits %>%
  select(trait_value_1.6.2, trait_value_1.6.3) %>%
  naniar::miss_var_summary(order = TRUE)
traits %>%
  select(trait_value_1.6.2, trait_value_1.6.3) %>%
  naniar::vis_miss(cluster = FALSE, sort_miss = TRUE)
```

``` r
traits %>%
  select(accepted_name, trait_value_1.6.2, trait_value_1.6.3) %>%
  filter(is.na(trait_value_1.6.2) | is.na(trait_value_1.6.3))
```

- Note the number of missing traits and calculate the ratio.

Check for synonyms

``` r
GIFT_species_lookup(genus = "Aconitum", epithet = "lycoctonum")
GIFT_species_lookup(genus = "Molinia", epithet = "arundinacea")
```

# 7 Save

``` r
write_csv(species, here("data", "processed", "data_processed_species.csv"))
write_csv(traits, here("data", "processed", "data_processed_traits.csv"))
```