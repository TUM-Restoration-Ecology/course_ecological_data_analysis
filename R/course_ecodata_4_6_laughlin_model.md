Analysis of Ecological Data <br> How to get functional traits?
================
<b>Markus Bauer</b> <br>
<b>2024-07-31</b>

- [Tasks](#tasks)
- [1 Preparation](#1-preparation)
  - [1.1 Load libraries and functions](#11-load-libraries-and-functions)
  - [1.2 Load data](#12-load-data)
- [2 Calculation of species
  compositions](#2-calculation-of-species-compositions)
  - [2.1 Using one trait](#21-using-one-trait)
    - [2.1.1 Constrain to CWM = 3.5](#211-constrain-to-cwm--35)
    - [2.1.2 Constrain to CWM = 3.5 and maximize Rao
      Q](#212-constrain-to-cwm--35-and-maximize-rao-q)
    - [2.1.3 Constrain to CWM = 3.5 and maximize Rao Q and
      Entropy](#213-constrain-to-cwm--35-and-maximize-rao-q-and-entropy)
    - [2.1.4 Maximize Rao Q and Entropy, no CWM trait
      constraint](#214-maximize-rao-q-and-entropy-no-cwm-trait-constraint)
  - [2.2 Using two traits](#22-using-two-traits)
    - [2.2.1 Constrain trait X to a CWM = 3.5, maximize Rao Q of trait
      Y](#221-constrain-trait-x-to-a-cwm--35-maximize-rao-q-of-trait-y)
    - [2.2.2 Constrain trait X to a CWM = 3.5, maximize Rao Q + Entropy
      of trait
      Y](#222-constrain-trait-x-to-a-cwm--35-maximize-rao-q--entropy-of-trait-y)
  - [2.3 Use the Ammer data](#23-use-the-ammer-data)
    - [2.3.1 Calculate species
      composition](#231-calculate-species-composition)
    - [2.3.2 Control species
      composition](#232-control-species-composition)
  - [2.4 Calculate abundances for several species
    compositions](#24-calculate-abundances-for-several-species-compositions)
    - [2.4.1 Taxonomic composition](#241-taxonomic-composition)
    - [2.4.2 Create a NMDS](#242-create-a-nmds)

**Markus Bauer**

Technichal University of Munich, TUM School of Life Sciences, Chair of
Restoration Ecology, Emil-Ramann-Straße 6, 85354 Freising, Germany

<markus1.bauer@tum.de>

<img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" alt="ORCID logo" width="20" height="20"/><https://orcid.org/0000-0001-5372-4174>

[![GoogleScholar](https://img.shields.io/badge/Google%20Scholar-4285F4?style=for-the-badge&logo=google-scholar&logoColor=white)](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/markus1bauer)

# Tasks

- Use in nmds
- Prepare Ammer data for analysis
- Select suitable values for SLA
- Use package Landscape metrics with CORINE landcover

[Original
vignette](https://cran.r-project.org/web/packages/Select/vignettes/selectSpecies.html)
(Laughlin et al. 2018)

# 1 Preparation

The package Select determines species probabilities (i.e., relative
abundances) that satisfy a given functional trait profile. Restoring
resilient ecosystems requires a flexible framework for selecting
assemblages that are based on the functional traits of species. However,
current trait-based models have been limited to algorithms that can only
select species by optimising specific trait values, and could not
elegantly accommodate the common desire among restoration ecologists to
produce functionally diverse assemblages. We have solved this problem by
applying a non-linear optimisation algorithm that optimises Rao Q, a
closed-form functional trait diversity index that incorporates species
abundances, subject to other linear constraints. This framework
generalises previous models that only optimised the entropy of the
community, and can optimise both functional diversity and entropy
simultaneously. This package can also be used to generate experimental
assemblages to test the effects of community-level traits on community
dynamics and ecosystem function.

## 1.1 Load libraries and functions

The Select package is of Laughlin et
al. ([2018](https://doi.org/10.1111/2041-210X.13023))

``` r
library(here)
```

    ## here() starts at C:/Users/marku/Documents/GitHub/course_ecological_data_analysis

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(Select)
rm(list = ls())
```

## 1.2 Load data

``` r
species_ammer <- read_csv(
  here("data", "processed", "data_processed_species.csv"),
  col_names = TRUE, col_types = cols(.default = "?")
  ) %>%
  column_to_rownames("accepted_name")
traits_ammer <- read_csv(
  here("data", "processed", "data_processed_traits.csv"),
  col_names = TRUE, col_types = cols(.default = "?")
  ) %>%
  column_to_rownames("accepted_name") %>%
  select(
    -name, -taxonomic_status, -starts_with("accepted_"), -starts_with("cv_"), -starts_with("n_")
    )
# sites_ammer <- read_csv(
#   here("data", "processed", "data_processed_sites.csv"),
#   col_names = TRUE, col_types = cols(.default = "?")
#   )
```

``` r
# herbs <- traits %>%
#   filter(!(family == "Poaceae" | family == "Fabaceae")) %>%
#   select(name, family)
# grass <- traits %>%
#   filter(family == "Poaceae") %>%
#   select(name, family)
# legumes <- traits %>%
#   filter(family == "Fabaceae") %>%
#   select(name, family)
```

# 2 Calculation of species compositions

## 2.1 Using one trait

Create a simple trait dataset using a species pool of 5 species with
trait values 1 through 5.

``` r
traits_one <- tibble(.rows = 5) %>%
  mutate(
    traitX = 1:5, # Set trait values for each species from 1 to 5
    name = letters[row_number()] # Most datasets have species names as row names, so let's add arbitrary row names using letters
    ) %>%
  column_to_rownames("name") %>%
  as.matrix() # Note that we make this a matrix to pass it into the function
```

### 2.1.1 Constrain to CWM = 3.5

Let us start with the most basic use of the function: to derive a
species abundance distribution where we only want to constrain the
abundances so that the community has a community-weighted mean (CWM)
trait equal to 3.5.

We will define four arguments in the function for this example t2c: this
is the matrix of species trait values that we want to constrain
constraints: this is a vector of CWM trait values. In this case, we only
have a one dimensional matrix of t2c, so this vector should contain only
one element: 3.5 t2d: when you are not maximizing functional diversity,
simply specify the same matrix here as you did for t2c obj: this is the
objective function that is being maximized. In this example we are not
maximizing functional diversity, so we use the entropy function (H)

``` r
result1 <- Select::selectSpecies(
  t2c = traits_one,
  constraints = c(traitX = 3.5),
  t2d = traits_one,
  obj = "H"
  )
```

    ## 
    ## Iter: 1 fn: -1.5461   Pars:  0.11205 0.14490 0.18738 0.24231 0.31335
    ## Iter: 2 fn: -1.5461   Pars:  0.11205 0.14490 0.18738 0.24231 0.31335
    ## solnp--> Completed in 2 iterations

``` r
Select::plotProbs(result1, traits_one, xlab = "Species")
```

![](course_ecodata_3_6_laughlin_model_files/figure-gfm/example1-1.png)<!-- -->

### 2.1.2 Constrain to CWM = 3.5 and maximize Rao Q

We will now see what happens when we maximize Rao Q (Q, quadratic
entropy), an index of functional diversity, rather than maximizing
entropy.

``` r
result2 <- Select::selectSpecies(
  t2c = traits_one,
  constraints = c(traitX = 3.5),
  t2d = traits_one,
  obj = "Q"
  )
# Note the only difference with result1 is a different objective function (Q).

Select::plotProbs(result2, traits_one, xlab = "Species")
```

Interestingly, the CWM trait value is the same for both result1 and
result2, but the species abundance distributions are radically
different. When maximizing Rao Q, this makes the most functionally
dissimilar species the most abundant, and all species in the middle of
the trait distribution have vanishingly small abundances. This is not a
desirable solution for ecological restoration. One way to fix this is to
optimize both Rao Q and Entropy simultaneously.

### 2.1.3 Constrain to CWM = 3.5 and maximize Rao Q and Entropy

We will now see what happens when we maximize a function that additively
combines both Rao Q (Q, quadratic entropy) and entropy (H).

``` r
result3 <- Select::selectSpecies(
  t2c = traits_one,
  constraints = c(traitX = 3.5),
  t2d = traits_one,
  obj="QH"
  )
### Note the only difference with result2 is a different objective function (QH)

Select::plotProbs(result3, traits_one, xlab = "Species")
```

Note that the abundance distribution still maximizes the most dissimilar
species, but it evens out the abundances across all the species.

### 2.1.4 Maximize Rao Q and Entropy, no CWM trait constraint

Suppose we do not want to constrain the abundances to satisfy a specific
CWM trait value, and simply want to maximize functional diversity. If
you do not want to constrain the results to satisfy a particular CWM
trait value, then leave the t2c and constraints arguments blank.

``` r
result4 <- Select::selectSpecies(
  t2d = traits_one,
  obj = "QH"
  )

Select::plotProbs(result4, traits_one, xlab = "Species")
```

## 2.2 Using two traits

In many cases, we will want to restore ecological communities with
convergence toward one trait value, but we want to diversify a different
trait. The following examples illustrate how to do so on a dataset with
known structure to easily illustrate the results.

Create 2-dimensional trait matrix

``` r
### 2-dimensional trait matrix of 16 species, evenly spaced between trait values 1 through 4
traits_two <- tibble(.rows = 16) %>%
  mutate(
    traitX = c(rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4)),
    traitY = c(rep(c(1, 2, 3, 4), 4)),
    name = letters[row_number()]
  ) %>%
  column_to_rownames("name")
traits_x <- traits_two %>%
  select(traitX) %>%
  as.matrix()
traits_y <- traits_two %>%
  select(traitY) %>%
  as.matrix()
```

### 2.2.1 Constrain trait X to a CWM = 3.5, maximize Rao Q of trait Y

``` r
result5 <- Select::selectSpecies(
  t2c = traits_x,
  constraints = c(traitX = 3.5),
  t2d = traits_y,
  obj = "Q",
  capd = FALSE
  )

Select::plotProbs(result5, traits_two, cex.lab = 0.7)
```

Note how this result is not desirable because it suppresses the
abundances of species with intermediate trait Y values. This is because
we only maximized Rao Q. Let us see how the results change when we
maximize both Rao Q and entropy.

### 2.2.2 Constrain trait X to a CWM = 3.5, maximize Rao Q + Entropy of trait Y

``` r
result6 <- Select::selectSpecies(
  t2c = traits_x,
  constraints = c(traitX = 3.5),
  t2d = traits_y,
  obj = "QH",
  capd = TRUE
  )

Select::plotProbs(result6, traits_two, cex.lab = 0.7)
```

## 2.3 Use the Ammer data

### 2.3.1 Calculate species composition

You can calculate with the Ammer data three different species
compositions.

First, calculate a species composition with a low CWM of SLA.

``` r
# plotcompData <- compData[which(compData$comp == i), ]
# row.names(plotcompData) <- plotcompData[, "name"]
# plotcompData <- plotcompData[, -(1:2)]
# mix_low <- Select::selectSpecies(
#   t2c = as.matrix(plotcompData),
#   constraints = c(sla = 2.995732, grass = 0.3, legume = 0.15),
#   t2d = seedmass,
#   obj = "QH",
#   capd = TRUE
# )
# ratioResults <- append(ratioResults, mix$prob)
# Select::plotProbs(mix_low, trait, xlab = "Species")
```

Second, calculate a species composition with a high CWM of SLA.

``` r
# plotcompData <- compData[which(compData$comp == i), ]
# row.names(plotcompData) <- plotcompData[, "name"]
# plotcompData <- plotcompData[, -(1:2)]
# mix_high <- Select::selectSpecies(
#   as.matrix(plotcompData),
#   constraints = c(sla = 2.995732, grass = 0.3, legume = 0.15),
#   obj = "QH",
#   capd = TRUE
# )
# ratioResults <- append(ratioResults, mix$prob)
# Select::plotProbs(mix_high, trait, xlab = "Species")
```

Third, calculate a species composition with an intermediate CWM of SLA.

``` r
# plotcompData <- compData[which(compData$comp == i), ]
#   row.names(plotcompData) <- plotcompData[, "name"]
#   plotcompData <- plotcompData[, -(1:2)]
# mix_intermediate <- Select::selectSpecies(
#   as.matrix(plotcompData),
#   constraints = c(sla = 2.995732, grass = 0.3, legume = 0.15),
#   obj = "QH",
#   capd = TRUE
# )
# ratioResults <- append(ratioResults, mix$prob)
# Select::plotProbs(mix_intermediate, trait, xlab = "Species")
```

### 2.3.2 Control species composition

``` r
data <- result5$prob %>%
  as.data.frame() %>%
  rownames_to_column("name") %>%
  as_tibble() %>%
  mutate(across(where(is.numeric), ~ round(.x, digits = 3))) %>%
  arrange(desc(V1))
```

## 2.4 Calculate abundances for several species compositions

Sometimes, it is necessary to constrain the species pool to a certain
amount of species, which are a substitute of a larger species pool. It
can be necessary to replicate a functional composition by different
species compositions which are randomly town from the larger species
pool.

### 2.4.1 Taxonomic composition

``` r
# compData <- as.data.frame(replicate(15, {
#   comp <- c(sample(herbs$name, 12),
#             sample(grass$name, 5),
#             sample(legume$name, 3)
#             )
#   }
#   )) %>%
#   gather(compData, "comp", "name", 1:15)
# 
# table(compData$name)
# length(table(compData$name))
# 
# compData <- traits %>%
#   inner_join(compData, by = "name")
```

### 2.4.2 Create a NMDS

``` r
data <- result5$prob %>%
  as.data.frame()
```

Create own species compositions from your data

<a href="#top">Back to top</a>
