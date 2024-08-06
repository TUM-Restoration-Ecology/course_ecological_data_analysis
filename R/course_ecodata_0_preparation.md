Analysis of Ecological Data <br> Install R, RStudio, GitHub and GitHub
Desktop
================
<b>Markus Bauer & Christina Hartung</b> <br>
<b>2024-08-06</b>

- [Tasks](#tasks)
- [1 Why R?](#1-why-r)
- [2 Use eduroam](#2-use-eduroam)
- [3 Install R and RStudio](#3-install-r-and-rstudio)
  - [3.1 R](#31-r)
  - [3.2 RStudio Desktop](#32-rstudio-desktop)
- [4 Use RStudio](#4-use-rstudio)
  - [4.1 The tidyverse](#41-the-tidyverse)
  - [4.2 Run a test script](#42-run-a-test-script)
  - [4.3 Use the \`renv\`\` package](#43-use-the-renv-package)
- [5 Use GitHub](#5-use-github)
- [6 Install GitHub Desktop](#6-install-github-desktop)

**Markus Bauer**<sup>1</sup> & **Christina Hartung**<sup>2</sup>

<sup>1</sup> Technichal University of Munich, TUM School of Life
Sciences, Chair of Restoration Ecology, Emil-Ramann-Straße 6, 85354
Freising, Germany

<markus1.bauer@tum.de>

<img src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" alt="ORCID logo" width="20" height="20"/><https://orcid.org/0000-0001-5372-4174>

[![Google
Scholar](https://img.shields.io/badge/Google%20Scholar-4285F4?style=for-the-badge&logo=google-scholar&logoColor=white)](https://scholar.google.de/citations?user=oHhmOkkAAAAJ&hl=de&oi=ao)
[![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)](https://github.com/markus1bauer)

<sup>2</sup> University of Applied Sciences Weihenstephan-Triesdorf,
Institute of Ecology and Landscape, Am Hofgarten 1, Building A10, 85354
Freising, Germany

# Tasks

- 

# 1 Why R?

([John
Verzani](https://www.karlin.mff.cuni.cz/~kulich/vyuka/Rdoc/Verzani-SimpleR.pdf))

R is an open-source (GPL) statistical environment modeled after S and
S-Plus (<http://www.insightful.com>). The S language was developed in
the late 1980s at AT&T labs. The R project was started by Robert
Gentleman and Ross Ihaka of the Statistics Department of the University
of Auckland in 1995. It has quickly gained a widespread audience. It is
currently maintained by the R core-development team, a hard-working,
international team of *volunteer* developers. The R project web page
(<http://www.r-project.org>) is the main site for information on R. At
this site are directions for obtaining the software, accompanying
packages and other sources of documentation.

The benefits of R for an introductory student are

- R is free. R is open-source and runs on UNIX, Windows and Macintosh.
- R has an excellent built-in help system.
- R has excellent graphing capabilities.
- R’s language has a powerful, easy to learn syntax with many built-in
  statistical functions.
- The language is easy to extend with user-written functions.
- R is a computer programming language. For programmers it will feel
  more familiar than others and for new computer users, the next leap to
  programming will not be so large.

What is R lacking compared to other software solutions?

- It has a limited graphical interface (S-Plus has a good one). This
  means, it can be harder to learn at the outset.
- The command language is a programming language so students must learn
  to appreciate syntax issues etc.

# 2 Use eduroam

First, you should get access to the internet with eduroam. Here are
further information about eduroam ([TUM
FAQ](https://www.it.tum.de/it/eduroam/)) and how to get access
([BayernCollab](https://collab.dvb.bayern/display/TUMdocs/eduroam+Anleitungen)).

# 3 Install R and RStudio

## 3.1 R

The software R is maintained by an international team of developers who
make the language available through the web page of The Comprehensive R
Archive Network. The top of the web page provides three links for
downloading R. Follow the link that describes your operating system:
Windows, Mac, or Linux.

1.  Go to the download page of RStudio Desktop:
    <https://posit.co/download/rstudio-desktop/>
2.  Select the button ‘Download and install R’
3.  You are linked to the site of the Comprehensive R Archive Network:
    <https://cran.r-project.org>
4.  In the first box ‘Download and Install R’, click the Download button
    for your system. For example: ‘Download R for Windows’ or enter
    <https://cran.r-project.org/bin/windows/>
5.  Select the installation file for your system. Windows users should
    select the ‘base’ distribution.
6.  Run the installation file.

<a href="#top">Back to top</a>

## 3.2 RStudio Desktop

RStudio is a graphical user interface for R. RStudio is an application
like Microsoft Word—except that instead of helping you write in English,
RStudio helps you write in R.

1.  Go to the download page of RStudio Desktop:
    <https://posit.co/download/rstudio-desktop/>
2.  Click the button ‘Download RStudio Desktop for Windows’ or scroll
    down for other OS
3.  Run the installation file.

# 4 Use RStudio

Now that you have both R and RStudio on your computer, you can begin
using R by opening the RStudio program.

When you start RStudio you see the following in the console (left side):

``` r
# R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"
# Copyright (C) 2024 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64
# 
# R ist freie Software und kommt OHNE JEGLICHE GARANTIE.
# Sie sind eingeladen, es unter bestimmten Bedingungen weiter zu verbreiten.
# Tippen Sie 'license()' or 'licence()' für Details dazu.
# 
# R ist ein Gemeinschaftsprojekt mit vielen Beitragenden.
# Tippen Sie 'contributors()' für mehr Information und 'citation()',
# um zu erfahren, wie R oder R packages in Publikationen zitiert werden können.
# 
# Tippen Sie 'demo()' für einige Demos, 'help()' für on-line Hilfe, oder
# 'help.start()' für eine HTML Browserschnittstelle zur Hilfe.
# Tippen Sie 'q()', um R zu verlassen.
```

In the first line is the installed version of R. For the whole citation:

``` r
citation()
```

    ## To cite R in publications use:
    ## 
    ##   R Core Team (2024). _R: A Language and Environment for Statistical
    ##   Computing_. R Foundation for Statistical Computing, Vienna, Austria.
    ##   <https://www.R-project.org/>.
    ## 
    ## Ein BibTeX-Eintrag für LaTeX-Benutzer ist
    ## 
    ##   @Manual{,
    ##     title = {R: A Language and Environment for Statistical Computing},
    ##     author = {{R Core Team}},
    ##     organization = {R Foundation for Statistical Computing},
    ##     address = {Vienna, Austria},
    ##     year = {2024},
    ##     url = {https://www.R-project.org/},
    ##   }
    ## 
    ## We have invested a lot of time and effort in creating R, please cite it
    ## when using it for data analysis. See also 'citation("pkgname")' for
    ## citing R packages.

- Note for your material-and-methods section the R version you use (not
  the RStudio version) and the whole citation.

## 4.1 The tidyverse

You’ll also need to install some R packages. An R package is a
collection of functions, data, and documentation that extends the
capabilities of base R. Using packages is key to the successful use of
R. The majority of the packages that you will learn in this book are
part of the so-called tidyverse. The packages in the tidyverse share a
common philosophy of data and R programming, and are designed to work
together naturally.

You can install the complete tidyverse with a single line of code:

``` r
# install.packages("tidyverse")
```

On your computer, type that line of code in the console, and then press
enter to run it. R will download the packages from CRAN and install them
on your computer.

You will not be able to use the functions, objects, or help files in a
package until you load it with `library()`. Once you have installed a
package, you can load it using the `library()` function:

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

## 4.2 Run a test script

## 4.3 Use the \`renv\`\` package

`install.packages(c("tidyverse", "here", "renv"))`

# 5 Use GitHub

# 6 Install GitHub Desktop
