Analysis of Ecological Data <br> How to get functional traits?
================
<b>Markus Bauer</b> <br>
<b>2024-08-05</b>

- [1 Preparation](#1-preparation)
  - [1.1 Load libraries](#11-load-libraries)
  - [1.2 The example data set](#12-the-example-data-set)
  - [1.3 Load data](#13-load-data)
- [2 Calculate traits from measured
  values](#2-calculate-traits-from-measured-values)
- [3 Traits from databases](#3-traits-from-databases)
  - [3.1 Metadata of traits](#31-metadata-of-traits)
  - [3.2 Get trait values](#32-get-trait-values)
    - [3.2.1 Raw data and aggregated
      data](#321-raw-data-and-aggregated-data)
    - [3.2.2 Get the references of the
      measurements](#322-get-the-references-of-the-measurements)
    - [3.2.3 Download traits](#323-download-traits)
- [4 Name resolving](#4-name-resolving)
- [5 Merge](#5-merge)
- [6 Save](#6-save)
- [7 Your own dataset](#7-your-own-dataset)

**Markus Bauer**

Technichal University of Munich, TUM School of Life Sciences, Chair of
Restoration Ecology, Emil-Ramann-Straße 6, 85354 Freising, Germany

<markus1.bauer@tum.de>

<img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" alt="ORCID logo" width="20" height="20"/><https://orcid.org/0000-0001-5372-4174>

[![GoogleScholar](https://img.shields.io/badge/Google%20Scholar-4285F4?style=for-the-badge&logo=google-scholar&logoColor=white)](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/markus1bauer)

# 1 Preparation

## 1.1 Load libraries

You can always install missing libraries by `install.packages()` like
`install.packages("here")`. We use the `tidyverse` package (Wickham et
al. [2019](https://doi.org/10.21105/joss.01686)), the `GIFT` package
(Denelle et al. [2023](https://doi.org/10.1111/2041-210X.14213)) and the
package `TNRS` (Boyle et
al.[2013](https://doi.org/10.1186/1471-2105-14-16))

``` r
library(here)
library(tidyverse)
library(GIFT)
library(TNRS)
rm(list = ls())
```

Example how to get citation information

``` r
citation("here")
```

    ## Um Paket 'here' in Publikationen zu zitieren, nutzen Sie bitte:
    ## 
    ##   Müller K (2020). _here: A Simpler Way to Find Your Files_. R package
    ##   version 1.0.1, <https://CRAN.R-project.org/package=here>.
    ## 
    ## Ein BibTeX-Eintrag für LaTeX-Benutzer ist
    ## 
    ##   @Manual{,
    ##     title = {here: A Simpler Way to Find Your Files},
    ##     author = {Kirill Müller},
    ##     year = {2020},
    ##     note = {R package version 1.0.1},
    ##     url = {https://CRAN.R-project.org/package=here},
    ##   }

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

The dataset is also available on GitHub (Bauer et al. 2018)

[![DOI:10.14471/2018.38.006](http://img.shields.io/badge/DOI-10.14471/2018.38.006-informational.svg)](https://doi.org/10.14471/2018.38.006)

<https://github.com/markus1bauer/2018_alluvial_forest_river_ammer/tree/main>

## 1.3 Load data

``` r
species_ammer <- read_csv(
  here("data", "raw", "example_course_alluvial_forest_mb", "data_raw_species.csv"),
  col_names = TRUE, col_types = cols(.default = "?")
  )
```

<a href="#top">Back to top</a>

# 2 Calculate traits from measured values

A short example how to calculate functional traits by your own. For
example the specific leaf area (SLA) of trees (*Acer platanoides* and
*Tilia cordata*) in a greenhouse experiment in Dürnast (Bauer et
al. [2023](https://doi.org/10.1007/s00468-023-02391-8)).

Per tree, three leaves were chosen and the leaf area was measured and
latter the dry weight (Perez-Harguindeguy et
al. [2013](https://www.uv.es/jgpausas/papers/PerezHarguindeguy-2013-AJB_traits-handbook2.pdf),
section 3.1).

``` r
traits_bricks <- read_csv(
  here("data", "raw", "example_course_brick_trees_mb", "data_raw.csv"),
  col_names = TRUE, col_types = cols(.default = "?")
)
traits_bricks %>%
  select(plot, block, species, starts_with("leaf")) %>%
  head()
```

    ## # A tibble: 6 × 9
    ##   plot      block species leaf1Mass leaf2Mass leaf3Mass leaf1Area leaf2Area
    ##   <chr>     <chr> <chr>       <dbl>     <dbl>     <dbl>     <dbl>     <dbl>
    ## 1 AY30Amyac A     Acer         0.69      0.47      0.64     107.       82.5
    ## 2 AY30Anoac A     Acer         0.24      0.3       0.26      40.2      53.3
    ## 3 AY30Tmyac A     Tilia        0.18      0.15      0.15      34.4      32.3
    ## 4 AY30Tnoac A     Tilia        0.12      0.11      0.09      24.4      22.3
    ## 5 BY30Amyac B     Acer         0.42      0.75      0.61      65.1      92.7
    ## 6 BY30Anoac B     Acer         0.71      0.46      0.57     102.       73.4
    ## # ℹ 1 more variable: leaf3Area <dbl>

Divide the leaf area by the leaf mass:

``` r
data <- traits_bricks %>%
  mutate(
    sla1 = leaf1Area / leaf1Mass,
    sla2 = leaf2Area / leaf2Mass,
    sla3 = leaf3Area / leaf3Mass,
    sla_mean = (sla1 + sla2 + sla3) / 3
  ) %>%
  select(1:3, 8, starts_with("leaf"), starts_with("sla"))
data %>%
  select(plot, species, starts_with("sla"))
```

    ## # A tibble: 100 × 6
    ##    plot      species  sla1  sla2  sla3 sla_mean
    ##    <chr>     <chr>   <dbl> <dbl> <dbl>    <dbl>
    ##  1 AY30Amyac Acer     155.  176.  176.     169.
    ##  2 AY30Anoac Acer     168.  178.  169.     171.
    ##  3 AY30Tmyac Tilia    191.  215.  193.     200.
    ##  4 AY30Tnoac Tilia    203.  203.  244.     217.
    ##  5 BY30Amyac Acer     155   124.  143.     141.
    ##  6 BY30Anoac Acer     144.  160.  147.     150.
    ##  7 BY30Tmyac Tilia    186.  164.  171      173.
    ##  8 BY30Tnoac Tilia    177.  161.  190      176.
    ##  9 CY30Amyac Acer     147.  152.  144.     148.
    ## 10 CY30Anoac Acer     141.  149.  141.     143.
    ## # ℹ 90 more rows

<a href="#top">Back to top</a>

# 3 Traits from databases

[Original
vignette](https://biogeomacro.github.io/GIFT/articles/GIFT.html#trait-data)
(Denelle & Weigelt 2024)

Get the functional trait values from the GLobal Inventory of Floras and
Traits (GIFT) database (Weigelt et
al. [2023](https://doi.org/10.1111/jbi.13623)) with the package `GIFT`
(Denelle et al. [2023](https://doi.org/10.1111/2041-210X.14213)).

## 3.1 Metadata of traits

There are many functional traits available in the GIFT database. Each of
these traits has an identification number called `trait_ID`. Since the
two functions for retrieving trait values, `GIFT_traits()` and
`GIFT_traits_raw()`, rely on these IDs, the first step is to call the
function `GIFT_traits_meta()` to know what the ID of the desired trait
is. For example, let’s say we want to retrieve the maximum vegetative
heights of plant species.

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

<a href="#top">Back to top</a>

## 3.2 Get trait values

### 3.2.1 Raw data and aggregated data

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
height <- GIFT::GIFT_traits(
  trait_IDs = c("1.6.2"), agreement = 0.66, bias_ref = FALSE, bias_deriv = FALSE
  )

height_raw <- GIFT::GIFT_traits_raw(trait_IDs = c("1.6.2"))

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

<a href="#top">Back to top</a>

### 3.2.2 Get the references of the measurements

If you want to look up the references that led to the aggregated trait
value, you can run this chunk:

``` r
references <- GIFT::GIFT_references(GIFT_version = "latest")

unique(unlist(strsplit(height$references_1.6.2, ",")))

references <- references[
  which(
    references$ref_ID %in% unique(
      unlist(strsplit(height$references_1.6.2, ","))
      )
    ), ]
references[1:2, ]
```

<a href="#top">Back to top</a>

### 3.2.3 Download traits

Here is an example how to download several traits:

``` r
trait_ids <- c(
  "1.2.2", "1.5.1", "1.6.3", "3.2.3", "3.3.1", "3.6.1", "3.7.1", "3.21.1",
  "4.1.3", "4.5.1"
  )

trait_meta %>%
  filter(Lvl3 %in% trait_ids) # Get an overview of selected traits

data_gift <- GIFT::GIFT_traits(
  trait_IDs = trait_ids,
  agreement = 0.66, bias_ref = FALSE, bias_deriv = FALSE
  )

rm(list = setdiff(ls(), c("species_ammer", "traits_ammer", "data_gift")))
```

<a href="#top">Back to top</a>

# 4 Name resolving

[Original
vignette](https://cran.r-project.org/web/packages/TNRS/vignettes/TNRS_vignette.html)
(Maitner 2024)

Use the package TNRS (Boyle et
al.[2013](https://doi.org/10.1186/1471-2105-14-16)) for name resolving
and the taxonoic resouces World Checklist of Vascular Plants of Kew
Gardens (Govaerts
[2023](http://sftp.kew.org/pub/data-repositories/WCVP/)) and World Flora
Online (WFO Consortium [2023](https://doi.org/10.5281/zenodo.8079052))

``` r
metadata <- TNRS::TNRS_metadata()

data <- species_ammer %>%
  rowid_to_column("id") %>%
  select(id, name) %>%
  TNRS::TNRS(
    sources = c("wcvp", "wfo"),
    classification = "wfo", # family classification
    mode = "resolve"
  )
names <- data %>%
  select(
    Name_submitted, Taxonomic_status, Accepted_name, Accepted_name_url,
    Accepted_family
    ) %>%
  rename_with(tolower)
```

<a href="#top">Back to top</a>

# 5 Merge

Merge the resolved species names and the functional plant traits into
the existing tables ‘species’ and ‘traits’.

``` r
data <- species_ammer %>%
  rename(name_submitted = name) %>%
  full_join(
    names %>% select(name_submitted, accepted_name), by = "name_submitted"
    )
species_ammer <- data

data <- names %>%
  left_join(
    data_gift %>% rename(accepted_name = work_species), by = "accepted_name"
    )
traits_ammer <- data %>%
  select(-work_ID, -work_author, -starts_with("references"))

rm(list = setdiff(ls(), c("species_ammer", "traits_ammer")))
```

<a href="#top">Back to top</a>

# 6 Save

``` r
write_csv(
  species_ammer, here("data", "processed", "data_processed_species_ammer.csv")
  )
write_csv(
  traits_ammer, here("data", "processed", "data_processed_traits_ammer.csv")
  )
```

<a href="#top">Back to top</a>

# 7 Your own dataset

- Select at least five functional plant traits.

- Note the categories, IDs, unit and the version of the GIFT database
  for the material-and-methods section.

- Get aggregated trait values from the GIFT database

- Resolve names with TNRS

- Note the version of the database from `metadata$version` and the
  sources of the names with `metadata$citations` for the
  material-and-methods section

- Note the range of ‘Overall_score’ for your submitted species names.
