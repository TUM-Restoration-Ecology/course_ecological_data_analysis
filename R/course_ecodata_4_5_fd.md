Analysis of Ecological Data <br> 3.5 Functional Diversity (FD)
================
<b>Markus Bauer</b> <br>
<b>2024-07-31</b>

- [Tasks](#tasks)
- [1 Preparation](#1-preparation)
  - [1.1 Load libraries and functions](#11-load-libraries-and-functions)
  - [1.2 Load data](#12-load-data)
- [2 Calculation of FD](#2-calculation-of-fd)
  - [2.1 Calculation with the dbFD
    function](#21-calculation-with-the-dbfd-function)
  - [2.2 Spain data](#22-spain-data)
    - [2.2.1 Numeric traits](#221-numeric-traits)
    - [2.2.2 FD with fuzzy coding](#222-fd-with-fuzzy-coding)
  - [2.3 Ammer data](#23-ammer-data)
- [3 Calculation of Alpha, beta and gamma
  FD](#3-calculation-of-alpha-beta-and-gamma-fd)

**Markus Bauer**

Technichal University of Munich, TUM School of Life Sciences, Chair of
Restoration Ecology, Emil-Ramann-Straße 6, 85354 Freising, Germany

<markus1.bauer@tum.de>

<img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" alt="ORCID logo" width="20" height="20"/><https://orcid.org/0000-0001-5372-4174>

[![Google
Scholar](https://img.shields.io/badge/Google%20Scholar-4285F4?style=for-the-badge&logo=google-scholar&logoColor=white)](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/markus1bauer)

# Tasks

- Make script tidy
- Quellen verlinken
- Show example Ammer
- Plot results of traits data
- Include Zeleny

This script uses large parts from [Bello et
al. (2021)](https://doi.org/10.1017/9781108628426)

# 1 Preparation

In this exercise we will learn how to compute different indices of
Functional diversity (FD). This exercise follows the second part of
Chapter 5 of the reference textbook, so all theoretical and mathematical
issues beyond the indices described in this exercise can be found there.
We will work with invented data and also field data from an climatic
gradient in NE Spain. In Chapter 3, and the relative exercise, we
already explained how to compute trait dissimilarity between species
pairs, and thus in this exercise we assume that users have already an
idea how to compute and interpret such trait dissimilarity. This file
will first cover the use of the `dbFD` function, with invented data and
then with the NE Spain data. Finally, we will learn how to compute
alpha, beta and gamma functional diversity with the Rao index.

## 1.1 Load libraries and functions

Load the libraries.

``` r
library(here)
library(tidyverse)
library(FD)
rm(list = ls())
```

We will also need one ad-hoc function ‘Rao’

``` r
source(here::here("data", "raw", "bello_etal-2021",  "chapter5", "Rao.r"))
```

## 1.2 Load data

``` r
species_spain <- read_delim(
  here::here("data", "raw", "bello_etal-2021", "chapter5", "speciesXplotsNE.txt"),
  col_names = TRUE, delim = "\t", col_types = cols(.default = "?")
  ) %>%
  column_to_rownames("species")
traits_spain <- read_delim(
  here::here("data", "raw", "bello_etal-2021", "chapter5", "speciesXtraitsNE.txt"),
  col_names = TRUE, delim = "\t", col_types = cols(.default = "?")
  ) %>%
  column_to_rownames("species")

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
```

# 2 Calculation of FD

## 2.1 Calculation with the dbFD function

As well summarized in the help function of `?dbFD`, dbFD implements a
flexible distance-based framework to compute multidimensional functional
diversity (FD) indices. dbFD returns the three FD indices of Villéger et
al. ([2008](https://doi.org/10.1890/07-1206.1)): functional richness
(FRic), functional evenness (FEve), and functional divergence (FDiv), as
well as functional dispersion (FDis; Laliberté and Legendre
[2010](https://doi.org/10.1890/08-2244.1)), Rao’s quadratic entropy (Q)
(Botta-Dukát
[2005](https://doi.org/10.1111/j.1654-1103.2005.tb02393.x)), a
functional group richness (FGR) (number of functional groups, see
Chapter 3 and Reference material Ch3 and Petchey and Gaston,
[2006](https://doi.org/10.1111/j.1461-0248.2006.00924.x)), and the
community-level weighted means of trait values (CWM; e.g. Lavorel et
al. [2008](https://doi.org/10.1111/j.1365-2435.2007.01339.x)). Some of
these FD indices consider species abundances (see Fig. 5.4 in the
reference book).

As we just saw, the function `dbFD` is very practical because it
computes a lot of indices of FD at once. Most of these indices are
introduced and discussed in the reference textbook, in Chapter 5. We
will thus learn how to apply the function `dbFD`, using first some of
the examples available in the help of the function, i.e. from `?dbFD`.
As in the exercises with CWM, the examples for the function consider
invented data, the ‘dummy’ data, containing both a ‘species x community’
matrix and a ‘species x trait’ matrix.

``` r
dummy$abun # Species x community matrix
```

    ##       sp1 sp2 sp3 sp4 sp5 sp6 sp7 sp8
    ## com1    1   1   0   0   4   2   0   0
    ## com2    0   0   0   2   1   0   0   5
    ## com3    2   0   0   0   0   1   0   3
    ## com4    1   0   7   0   0   0   0   0
    ## com5    0   0   2   3   3   0   0   0
    ## com6    0   3   0   0   5   6   1   6
    ## com7    3   5   0   3   0   0   0   0
    ## com8    0   0   0   0   6   2   1   2
    ## com9    4   1   1   3   0   0   2   0
    ## com10   0   4   1   0   0   0   6   1

``` r
dummy$trait # Species x trait matrix
```

    ##     num1 num2 fac1 fac2 ord1 ord2 bin1 bin2
    ## sp1  9.0  4.5    A    X    3    2    0    1
    ## sp2  8.1  6.0    A    Z <NA>    1    0    1
    ## sp3   NA  2.3    C    Y    5    3    1    1
    ## sp4  3.2  5.4    B    Z    1    7    0    0
    ## sp5  5.8  1.2    C    X    2    6   NA    0
    ## sp6  3.4  8.5    C    Y    2    1    1    1
    ## sp7  7.5  2.1    B    X    3    2    1    0
    ## sp8  4.3  6.5 <NA>    Z    1    3    0    1

*Notice* that the function `dbFD`, while extremely useful, has also many
strict requirements. The species x community matrix should be a data
frame and the order and name of species should be *exactly* the same in
the ‘species x community matrix’ and in the Species x trait matrix.

It is very important to check that the species that you have in these
‘species x community’ matrix and ‘species x trait’ matrices should be
the same. In this case we know that they are the same, so we can already
run the code:

``` r
ex1 <- dbFD(dummy$trait, dummy$abun)
```

    ## Species x species distance matrix was not Euclidean. 'sqrt' correction was applied. 
    ## FEVe: Could not be calculated for communities with <3 functionally singular species. 
    ## FRic: To respect s > t, FRic could not be calculated for communities with <3 functionally singular species. 
    ## FRic: Dimensionality reduction was required. The last 5 PCoA axes (out of 7 in total) were removed. 
    ## FRic: Quality of the reduced-space representation (based on corrected distance matrix) = 0.6108352 
    ## FDiv: Could not be calculated for communities with <3 functionally singular species.

First we see that the function is very “talkative”, i.e. it explains a
lot of things while running. We can eventually switch-off the messages
by using `messages = F`, as we will do in the next examples. In the
meantime we can see that FEve (functional evenness) cannot be computed
in communities with less than 3 species. Similarly as FRic (functional
richness) and FDiv (functional divergence) cannot be computed if there
is less that 3 unique species, meaning that they should at least have
one different trait value. Other important messages include the number
of PCoA axes considered. Out of the 7 considered (i.e. with 8 species,
the number of PCoA is 7) a total of 2 axes are retained to reflect trait
differences between species (see Chapter 3 in the reference book and R
material Ch 3. In other words the PCoA has synthesized all traits into
two multivariate axes.

Let’s now see the results. The results are stored in the object ex1,
which is a ‘list’. As such if you need a particular object you can
extract it, for example as `ex1$FRic`:

``` r
class(ex1)
```

    ## [1] "list"

``` r
ex1
```

    ## $nbsp
    ##  com1  com2  com3  com4  com5  com6  com7  com8  com9 com10 
    ##     4     3     3     2     3     5     3     4     5     4 
    ## 
    ## $sing.sp
    ##  com1  com2  com3  com4  com5  com6  com7  com8  com9 com10 
    ##     4     3     3     2     3     5     3     4     5     4 
    ## 
    ## $FRic
    ##        com1        com2        com3        com4        com5        com6 
    ## 0.174201349 0.102097174 0.002157642          NA 0.143151204 0.231083703 
    ##        com7        com8        com9       com10 
    ## 0.073683375 0.174220613 0.337446205 0.228904118 
    ## 
    ## $qual.FRic
    ## [1] 0.6108352
    ## 
    ## $FEve
    ##      com1      com2      com3      com4      com5      com6      com7      com8 
    ## 0.8432334 0.4628635 0.8659657        NA 0.9081209 0.8651734 0.7681991 0.7606650 
    ##      com9     com10 
    ## 0.7994638 0.4944601 
    ## 
    ## $FDiv
    ##      com1      com2      com3      com4      com5      com6      com7      com8 
    ## 0.8429221 0.8593250 0.6303031        NA 0.7631346 0.8966687 0.8356095 0.8681163 
    ##      com9     com10 
    ## 0.7015118 0.9712554 
    ## 
    ## $FDis
    ##      com1      com2      com3      com4      com5      com6      com7      com8 
    ## 0.3481687 0.1670560 0.2375808 0.1146261 0.3211366 0.3302330 0.2532751 0.2877931 
    ##      com9     com10 
    ## 0.3421687 0.3503927 
    ## 
    ## $RaoQ
    ##       com1       com2       com3       com4       com5       com6       com7 
    ## 0.12835440 0.04063622 0.06497224 0.03003235 0.11858388 0.11738479 0.07868047 
    ##       com8       com9      com10 
    ## 0.09308505 0.12431265 0.12490783 
    ## 
    ## $CWM
    ##           num1     num2 fac1 fac2 ord1 ord2 bin1 bin2
    ## com1  5.887500 4.037500    C    X    2    6    0    0
    ## com2  4.212500 5.562500    B    Z    1    3    0    1
    ## com3  5.716667 6.166667    A    Z    1    3    0    1
    ## com4  9.000000 2.575000    C    Y    5    3    1    1
    ## com5  4.500000 3.050000    C    Z    1    7    0    0
    ## com6  5.095238 5.528571    C    Z    2    1    0    1
    ## com7  7.009091 5.427273    A    Z    1    1    0    1
    ## com8  5.245455 3.572727    C    X    2    6    1    0
    ## com9  6.870000 4.245455    A    X    3    2    0    1
    ## com10 7.427273 3.783333    B    X    3    2    1    0

``` r
ex1$FRic
```

    ##        com1        com2        com3        com4        com5        com6 
    ## 0.174201349 0.102097174 0.002157642          NA 0.143151204 0.231083703 
    ##        com7        com8        com9       com10 
    ## 0.073683375 0.174220613 0.337446205 0.228904118

For categorical and binary traits you might want to include the argument
`CWM.type = "all"` when running it.

``` r
ex1 <- dbFD(dummy$trait, dummy$abun, CWM.type = "all", message = F)
# notice we add `message = F` to avoid too many annoying messages from the function :)
ex1$CWM
```

    ##           num1     num2    fac1_A     fac1_B     fac1_C    fac2_X     fac2_Y
    ## com1  5.887500 4.037500 0.2500000 0.00000000 0.75000000 0.6250000 0.25000000
    ## com2  4.212500 5.562500 0.0000000 0.25000000 0.12500000 0.1250000 0.00000000
    ## com3  5.716667 6.166667 0.3333333 0.00000000 0.16666667 0.3333333 0.16666667
    ## com4  9.000000 2.575000 0.1250000 0.00000000 0.87500000 0.1250000 0.87500000
    ## com5  4.500000 3.050000 0.0000000 0.37500000 0.62500000 0.3750000 0.25000000
    ## com6  5.095238 5.528571 0.1428571 0.04761905 0.52380952 0.2857143 0.28571429
    ## com7  7.009091 5.427273 0.7272727 0.27272727 0.00000000 0.2727273 0.00000000
    ## com8  5.245455 3.572727 0.0000000 0.09090909 0.72727273 0.6363636 0.18181818
    ## com9  6.870000 4.245455 0.4545455 0.45454545 0.09090909 0.5454545 0.09090909
    ## com10 7.427273 3.783333 0.3333333 0.50000000 0.08333333 0.5000000 0.08333333
    ##          fac2_Z     ord1_1    ord1_2     ord1_3     ord1_5     ord2_1
    ## com1  0.1250000 0.00000000 0.7500000 0.12500000 0.00000000 0.37500000
    ## com2  0.8750000 0.87500000 0.1250000 0.00000000 0.00000000 0.00000000
    ## com3  0.5000000 0.50000000 0.1666667 0.33333333 0.00000000 0.16666667
    ## com4  0.0000000 0.00000000 0.0000000 0.12500000 0.87500000 0.00000000
    ## com5  0.3750000 0.37500000 0.3750000 0.00000000 0.25000000 0.00000000
    ## com6  0.4285714 0.28571429 0.5238095 0.04761905 0.00000000 0.42857143
    ## com7  0.7272727 0.27272727 0.0000000 0.27272727 0.00000000 0.45454545
    ## com8  0.1818182 0.18181818 0.7272727 0.09090909 0.00000000 0.18181818
    ## com9  0.3636364 0.27272727 0.0000000 0.54545455 0.09090909 0.09090909
    ## com10 0.4166667 0.08333333 0.0000000 0.50000000 0.08333333 0.33333333
    ##           ord2_2     ord2_3    ord2_6    ord2_7    bin1_0    bin1_1    bin2_0
    ## com1  0.12500000 0.00000000 0.5000000 0.0000000 0.2500000 0.2500000 0.5000000
    ## com2  0.00000000 0.62500000 0.1250000 0.2500000 0.8750000 0.0000000 0.3750000
    ## com3  0.33333333 0.50000000 0.0000000 0.0000000 0.8333333 0.1666667 0.0000000
    ## com4  0.12500000 0.87500000 0.0000000 0.0000000 0.1250000 0.8750000 0.0000000
    ## com5  0.00000000 0.25000000 0.3750000 0.3750000 0.3750000 0.2500000 0.7500000
    ## com6  0.04761905 0.28571429 0.2380952 0.0000000 0.4285714 0.3333333 0.2857143
    ## com7  0.27272727 0.00000000 0.0000000 0.2727273 1.0000000 0.0000000 0.2727273
    ## com8  0.09090909 0.18181818 0.5454545 0.0000000 0.1818182 0.2727273 0.6363636
    ## com9  0.54545455 0.09090909 0.0000000 0.2727273 0.7272727 0.2727273 0.4545455
    ## com10 0.50000000 0.16666667 0.0000000 0.0000000 0.4166667 0.5833333 0.5000000
    ##          bin2_1
    ## com1  0.5000000
    ## com2  0.6250000
    ## com3  1.0000000
    ## com4  1.0000000
    ## com5  0.2500000
    ## com6  0.7142857
    ## com7  0.7272727
    ## com8  0.3636364
    ## com9  0.5454545
    ## com10 0.5000000

In the second example of the help page for `dbFD`, it is possible to see
a the argument `w`, which was discussed for the `gowdis` in the R
material Ch3 and in Chapter 3 of the reference textbook. As a matter of
fact the function gowdis is also included in the `dbFD` function and,
although you do not see it, it is applied exactly as already introduced
in the R material 3. This is the first step for the calculation of
functional diversity.

The argument `w` allows us to give different weights to the traits, for
example when you want some traits to be more important in the
calculation of the dissimilarity between species and thus functional
diversity (we discuss the importance of this weight in the exercises
related to Chapter 3). The example 2 in the help function looks like
this:

``` r
# add variable weights
w <- c(1, 5, 3, 2, 5, 2, 6, 1)
ex2 <- dbFD(dummy$trait, dummy$abun, w, corr = "cailliez", message = FALSE)
# `cailliez` correction is used because `sqrt` does not work
```

In practice the second trait got 5 times more weight than the first one
(the vector w, includes one value for each trait, reflecting the
intended weight in the calculation; again see R material Ch 2.

It might happen that when running the dbFD function you get some sort of
errors. For example if you run a line such as
`ex2 <- dbFD(dummy$trait, dummy$abun, w, message = F)`, you will get one
saying that the distance between species did not have Euclidean
properties after the square root correction, which is automatically
applied in the function. The PCoA analyses work better when the
dissimilarity matrix (differences between species in terms of traits; R
material Ch 2) has Euclidean properties and the dbFD function will not
work unless this condition is met. In the case that the function does
not work directly with your data, you can try other corrections with the
argument corr. The help function provides enough information on this
issue. In practice the effect of such correction is not very strong, at
least in our experience, especially when there are enough species.

It is very important to notice that the function `dbFD` works not only
by providing a ‘species x trait’ matrix like the `dummy$trait` in the
example above. You can also provide directly a distance matrix
calculated before hand, for example using the Gower distance (see
Chapter 3; R material Ch 3). If you run the following example (ex3) you
will thus get the same results as in the object ex1.

``` r
trait.d <- gowdis(dummy$trait) # Gower distance
ex3 <- dbFD(trait.d, dummy$abun, message = FALSE)
ex1$FRic == ex3$FRic
```

    ##  com1  com2  com3  com4  com5  com6  com7  com8  com9 com10 
    ##  TRUE  TRUE  TRUE    NA  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE

Similarly if you apply `w` in the calculation of the Gower distance:

``` r
trait.dw <- gowdis(dummy$trait, w) # Gower distance with a different weight for the traits
ex3w <- dbFD(trait.dw, dummy$abun, corr = "cailliez", message = FALSE)
ex2$FRic == ex3w$FRic
```

    ##  com1  com2  com3  com4  com5  com6  com7  com8  com9 com10 
    ##  TRUE  TRUE  TRUE    NA  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE

Although we already discussed what is a dissimilarity matrix in Chapter
3, and its corresponding R material Ch 3, in the next lines, we show
once again how it does look:

``` r
trait.d # Gower distance
```

    ##           sp1       sp2       sp3       sp4       sp5       sp6       sp7
    ## sp2 0.2181884                                                            
    ## sp3 0.5240052 0.6678082                                                  
    ## sp4 0.6737443 0.5610028 0.8225701                                        
    ## sp5 0.5291113 0.8145699 0.4862253 0.4843264                              
    ## sp6 0.6100161 0.5932587 0.2784736 0.7073925 0.6067323                    
    ## sp7 0.4484235 0.6863374 0.4848663 0.5575126 0.3023416 0.6187844          
    ## sp8 0.4072834 0.2039443 0.5958904 0.2390962 0.5585525 0.4470207 0.7030186

In practice, this object shows the functional distance,
i.e. dissimilarity between each pair of species, in this case expressed
as an average dissimilarity over all (standardized) traits.

The function `dbFD` applies the Gower distance automatically if you use
a species x trait matrix (such as ‘spxt’). It will also work when
proving a dissimilarity matrix, unless such dissimilarity does not
contain any NAs (NAs in the dissimilarity matrix are not accepted!).
This means that you can either provide a specific dissimilarity matrix
after checking that there are no NAs or you can make sure that, in the
‘species x trait’ matrix, you have ‘enough’ trait values. In other
words, this means that species pairs should have, at least, information
for one common trait.

A special case occurs when only one trait is used, as a vector, and
there is some NA. In this case the function works (and it will remove
the species with missing values):

``` r
num1 <- dummy$trait[, 1] #take only one trait, as a vector, with missing values
names(num1) <- rownames(dummy$trait) #give "back" the species names
ex4 <- dbFD(num1, dummy$abun) #it works
```

    ## Warning: Species with missing trait values have been excluded. 
    ## FEVe: Could not be calculated for communities with <3 functionally singular species. 
    ## FDis: Equals 0 in communities with only one functionally singular species. 
    ## FRic: Only one continuous trait or dimension in 'x'. FRic was measured as the range, NOT as the convex hull volume. 
    ## FDiv: Cannot not be computed when 'x' contains one single continuous trait or dimension.

In all cases, you need to provide the species names in the trait data.
If you have only one trait, WITHOUT NAs, you can do

``` r
round(gowdis(dummy$trait["num2"]), 3) 
```

    ##       sp1   sp2   sp3   sp4   sp5   sp6   sp7
    ## sp2 0.205                                    
    ## sp3 0.301 0.507                              
    ## sp4 0.123 0.082 0.425                        
    ## sp5 0.452 0.658 0.151 0.575                  
    ## sp6 0.548 0.342 0.849 0.425 1.000            
    ## sp7 0.329 0.534 0.027 0.452 0.123 0.877      
    ## sp8 0.274 0.068 0.575 0.151 0.726 0.274 0.603

``` r
#does not remove species names and keep the trait data as a matrix; #here we keep only 3 decimals for visual purposes
ex.1trait.noNA <- dbFD(
  gowdis(dummy$trait[, "num2", drop = FALSE]), dummy$abun, message = FALSE
  ) #it works! no errors provided.
```

## 2.2 Spain data

### 2.2.1 Numeric traits

Let’s now apply the `dbFD` function to the NE Spain data already
described in section 5.1 above. The data was already loaded already
above. As in the exercise on CWM, we need to improve the normality of
the trait SLA in the matrix.

``` r
traits_spain$SLA <- log(traits_spain$SLA) #improve the normality of the trait values
head(traits_spain)
```

    ##          GrowhtForm LEG      SLA LF_Th LF_G LF_H LF_hCh LF_wCh LF_NP LF_P
    ## Acercamp      shrub   0 2.753661     0    0  0.0    0.0      0     0    1
    ## Achimill       forb   0 2.681022     0    0  0.5    0.5      0     0    0
    ## Aegigeni      grass   0 2.721295     1    0  0.0    0.0      0     0    0
    ## Alchhybr       forb   0 2.944439     0    0  1.0    0.0      0     0    0
    ## Anemhepa       forb   0 2.557227     0    0  1.0    0.0      0     0    0
    ## Anthmont       forb   1 2.602690     0    0  1.0    0.0      0     0    0

Let’s directly use the `dbFD` function. We need first to recall that we
have a fuzzy coding data in the trait matrix (`spxt$LF_xx`) and for now,
for simplicity, let’s forget about this type of trait data. We can use
thus only the first 3 columns of the trait matrix for the time being
(but see below for using all traits). Now that we know how the `dbFD`
function works, i.e. it is used very similarly to the function
functcomp, it is easy to have a lot of indices of FD computed in only
one line! But for the time being we avoid computing again the CWM
values, with the argument `calc.CWM = FALSE` and we also ask that the
FRic will be standardized between 0 and 1, i.e. `stand.FRic = TRUE`.

``` r
resFD <- dbFD(
  traits_spain[, 1:3], log(t(species_spain) + 1),
  message = FALSE, calc.CWM = FALSE, stand.FRic = TRUE
  ) 
```

    ## Warning in is.euclid(x.dist): Zero distance(s)
    ## Warning in is.euclid(x.dist): Zero distance(s)

    ## Warning in is.euclid(x.dist2): Zero distance(s)
    ## Warning in is.euclid(x.dist2): Zero distance(s)

    ## Warning in is.euclid(x.dist): Zero distance(s)
    ## Warning in is.euclid(x.dist): Zero distance(s)
    ## Warning in is.euclid(x.dist): Zero distance(s)

We can now explore a bit the results. First let’s see how much the
different FD indices are correlated between them:

``` r
important.indices <- cbind(
  resFD$nbsp, resFD$FRic, resFD$FEve, resFD$FDiv, resFD$FDis, resFD$RaoQ
  )
colnames(important.indices) <- c("NumbSpecies", "FRic", "FEve", "FDiv", "FDis", "Rao")
pairs(important.indices, pch = 20)
```

![](course_ecodata_3_5_fd_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

We generally see that FRic is positively correlated to the number of
species, as expected. Of course FDis and Rao give very similar results
(they are actually the same index, with only a squaring difference, see
Pavoine & Bonsall (2011)). Slight deviations from the perfect linear
relationship depends ONLY on the fact that Rao is computed directly from
a Gower trait distance, while in FDis the Gower trait distance is first
transformed into a PCoA and then the PCoA axes are used to compute the
dissimilarities in the multivariate space. In some cases FEve is
correlated negatively to FRic, as we discussed in Chapter 5 of the
reference text book (see below for more details). We also see that Rao
and FDis are increasing when the ‘range’ of traits, i.e. FRic (which is,
as we discuss in Chapter 5 of the reference text book, is the size of
the Convex Hull when multiple traits are considered, and the range, when
only one trait is considered). We also see that both FDis and Rao are
not correlated to the number of species. Finally, we also see that FDiv
is quite correlated with FDis and Rao.

### 2.2.2 FD with fuzzy coding

We remind users (see Chapter 3 of the reference textbook, and R material
Ch 3) that for computing FD with all traits considered in the ‘spxt’
matrix, which include fuzzy coding (‘LF_xx’ labeled column), we need to
complicate a bit more the script. Otherwise the Gower distance used in
the dbFD function will ‘understand’ that each of the LF_xx column in the
matrix is a different trait. We certainly do not want that, even if it
would provide nicer results, as we now want the 4 traits to have the
same weight and scale! We can use the approach already introduced for
the R material Ch 3, i.e. compute the dissimilarity for each trait
separately and then average the dissimilarity across traits:

``` r
head(traits_spain)
```

    ##          GrowhtForm LEG      SLA LF_Th LF_G LF_H LF_hCh LF_wCh LF_NP LF_P
    ## Acercamp      shrub   0 2.753661     0    0  0.0    0.0      0     0    1
    ## Achimill       forb   0 2.681022     0    0  0.5    0.5      0     0    0
    ## Aegigeni      grass   0 2.721295     1    0  0.0    0.0      0     0    0
    ## Alchhybr       forb   0 2.944439     0    0  1.0    0.0      0     0    0
    ## Anemhepa       forb   0 2.557227     0    0  1.0    0.0      0     0    0
    ## Anthmont       forb   1 2.602690     0    0  1.0    0.0      0     0    0

``` r
all.dist <- (gowdis(traits_spain["GrowhtForm"]) + gowdis(traits_spain["SLA"]) + gowdis(traits_spain["LEG"]) +
     gowdis(traits_spain[, 4:10]) / max(gowdis(traits_spain[, 4:10]))) / 4
resFD.alltraits <- dbFD(
  all.dist, log(t(species_spain) + 1), message = FALSE, calc.CWM = FALSE, stand.FRic = TRUE
  )
```

    ## Warning in is.euclid(x.dist): Zero distance(s)
    ## Warning in is.euclid(x.dist): Zero distance(s)

    ## Warning in is.euclid(x.dist2): Zero distance(s)
    ## Warning in is.euclid(x.dist2): Zero distance(s)

    ## Warning in is.euclid(x.dist): Zero distance(s)
    ## Warning in is.euclid(x.dist): Zero distance(s)
    ## Warning in is.euclid(x.dist): Zero distance(s)

## 2.3 Ammer data

``` r
# resFD <- dbFD(
#   species_ammer, traits_ammer, message = FALSE, calc.CWM = FALSE, stand.FRic = TRUE
#   ) 
```

# 3 Calculation of Alpha, beta and gamma FD

In Chapter 5 we also discuss options to partition functional diversity,
in a way similar to species diversity, i.e. in alpha and beta diversity.
The literature about biodiversity partitioning is generally very
complex, often with polarized views (e.g., Jost 2007). The partitioning
of functional diversity has not been solved for many indices. The most
reliable one seems Rao, following the work by de Bello et al. (2010) and
Botta-Dukát (2018). Other options involve computing convex hull, similar
to the FRic index, and look at intersections of the convex hull between
communities. This is done in the package betapart
(<https://cran.r-project.org/web/packages/betapart/betapart.pdf>) and
using the functions functional.beta.multi and functional.beta.pair. We
will introduce similar indice in the R material Ch 6. Here we focus
instead mainly on Rao, following the reference textbook Chapter 5 and
based on the fact that this approach can also consider species
abundances.

There has been a debate, in Journal of Ecology, between Villeger &
Mouillot (Villéger and Mouillot 2008) and Hardy and Jost (2008), on how
Rao should be expressed, which is summarized for non mathematicians in
the work by de Bello et al. (2010) . This work proposes solutions that
solves the debate mentioned above and applies the corrections proposed
by Jost for species diversity indices such as Simpson. All these
concepts are resolved in the R function Rao (F. de Bello, Lavergne, et
al. 2010), which is at the moment the only available function solving
these issues. We can open the function (which actually can also account
for phylogentic distance between species) and run it with few examples.
Notice that in this case the ‘species x communities’ matrix has species
as rows. We can first play with the simple example provided in Chapter 5
of the reference textbook, which we recreate below:

``` r
source(here::here("data", "raw", "bello_etal-2021",  "chapter5", "Rao.r"))
comm.1 <- c(1, 1, 0, 0, 0)
comm.2 <- c(0, 0, 1, 1, 0)
comm.3 <- c(0, 0, 1, 1, 1)
coms <- cbind(comm.1, comm.2, comm.3)
rownames(coms) <- paste("sp", 1:5, sep = ".")
coms
```

    ##      comm.1 comm.2 comm.3
    ## sp.1      1      0      0
    ## sp.2      1      0      0
    ## sp.3      0      1      1
    ## sp.4      0      1      1
    ## sp.5      0      0      1

``` r
trait.reg1 <- c("grass", "grass", "forb", "forb", "grass")
trait.reg2 <- c("grass", "grass", "legume", "forb", "grass")
partRao.Reg1 <- Rao(
  coms, dfunc = gowdis(as.data.frame(trait.reg1)),
  dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL
  )
partRao.Reg2 <- Rao(
  coms, dfunc = gowdis(as.data.frame(trait.reg2)),
  dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL
  )
```

The `Rao` function needs a species x communities matrix with species as
rows and does not need that the dissimilarity matrix will include
species names (it assumes, maybe wrongly, users are careful enough to
have species in the same order in the trait matrix and in the species
composition matrix). We just run the calculations, using the Jost
correction (`Jost = TRUE`) to avoid having an underestimation of beta
diversity (see Chapter 5 in the reference textbook for more details). We
can now find the results used which we need to build a figure similar to
Fig.5.6 in the reference textbook. The alpha taxonomic diversity (which
is the mean Simpson per plot), the beta and gamma, can be obtained as:

``` r
partRao.Reg1$TD$Mean_Alpha #alpha
```

    ## [1] 2.25

``` r
partRao.Reg1$TD$Beta_add #beta
```

    ## [1] 2.25

``` r
partRao.Reg1$TD$Mean_Alpha + partRao.Reg1$TD$Beta_add == partRao.Reg1$TD$Gamma #alpha+beta=gamma
```

    ## [1] TRUE

Indeed alpha + beta = gamma. We are aware some authors do not like to
express beta diversity in additivity terms, so it has to be expressed in
proportional terms, for example.

``` r
partRao.Reg1$TD$Beta_add * 100 / partRao.Reg1$TD$Gamma #beta / gamma
```

    ## [1] 50

``` r
partRao.Reg1$TD$Beta_prop #which is the same as above
```

    ## [1] 50

You can apply the same procedure between each pair of plot. The
taxonomical beta diversity between each pair of plots is shown in the
following object. Notice that between a pair of plots the value 50%
implies a complete functional replacement/turnover (see textbook for
more details). In fact, the value obtained, i.e. ‘Beta_prop’, depend on
the number of plots considered (the maximum turnover is 50% for two
plots and increases with increased number of plots). As such ‘Beta_prop’
can be standardized, to be independent on the number of plots, with the
formula 14 in de Bello et al. (2010).

``` r
round(partRao.Reg1$TD$Pairwise_samples$Beta_prop) #value obtained without standardization
```

    ##        comm.1 comm.2
    ## comm.2     50       
    ## comm.3     50     10

``` r
round(partRao.Reg1$TD$Pairwise_samples$Beta_prop / (1 - 1 / 2)) #standardized value, where 2 is the number of plos considered in the test, i.e. `turnover` between plots. 
```

    ##        comm.1 comm.2
    ## comm.2    100       
    ## comm.3    100     20

We can now try to make a figure similar to Fig.5.6.

``` r
TD <- c(partRao.Reg1$TD$Mean_Alpha, partRao.Reg1$TD$Beta_add)
FD.reg1 <- c(partRao.Reg1$FD$Mean_Alpha, partRao.Reg1$FD$Beta_add)
FD.reg2 <- c(partRao.Reg2$FD$Mean_Alpha, partRao.Reg2$FD$Beta_add)
barplot(cbind(TD, FD.reg1, FD.reg2), beside = F, col = c(8, 0),
        ylim = c(0, 5), ylab = c("Rao Equivalent Numbers"))
points(2, 4, pch = 22, bg = "grey")
points(2, 3.5, pch = 22, bg = "white")
text(2, 4, "alpha", pos = 4)
text(2, 3.5, "beta", pos = 4)
```

![](course_ecodata_3_5_fd_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

In practice we can see that in the second region (or in general when
considering the trait ‘trait.reg2’) we have higher alpha but also higher
beta diversity (i.e. functional replacement). As a proportion, however,
the functional replacement is the same in both regions (i.e., ~42% of
the trait diversity is between plots, as compared to within plots):

``` r
partRao.Reg1$FD$Beta_prop
```

    ## [1] 40.57971

``` r
partRao.Reg2$FD$Beta_prop
```

    ## [1] 42.42424

Let’s now work with the NE Spain data:

``` r
partRao.NESpain <- Rao(
  log(species_spain + 1), all.dist, dphyl = NULL, weight = FALSE, Jost = TRUE, structure = NULL
  )
partRao.NESpain$TD$Beta_prop
```

    ## [1] 74.17581

``` r
partRao.NESpain$FD$Beta_prop
```

    ## [1] 6.737896

``` r
TD.NEspain.prop <- c(100 * partRao.NESpain$TD$Mean_Alpha / partRao.NESpain$TD$Gamma,
                     partRao.NESpain$TD$Beta_prop)
FD.NEspain.prop <- c(100*partRao.NESpain$FD$Mean_Alpha/partRao.NESpain$FD$Gamma, 
                     partRao.NESpain$FD$Beta_prop)
barplot(cbind(TD.NEspain.prop, FD.NEspain.prop), beside = F, col = c(8, 0),
        ylim = c(0, 100), ylab = c("Proportion of Alpha and Beta in Equivalent numbers"))
text(0.7, 60, "beta")
text(0.7, 15, "alpha")
```

![](course_ecodata_3_5_fd_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

We can see in this figure that, despite the great turnover in species
composition (beta TD, grey colour, left stack) we have a real low
functional turnover. In other words **most of the trait dissimilarity
between species is found within a plot and not across plots**, even if
there are such a marked environmental changes and even if species
composition changes are very strong (high beta TD).
