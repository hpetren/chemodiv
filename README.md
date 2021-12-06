
<!-- README.md is generated from README.Rmd. Please edit that file -->

# chemdiv

<!-- badges: start -->
<!-- badges: end -->

chemdiv is an in development R package to analyse phytochemical data.

The package can be used to calculate various types of diversities and
dissimilarities for phytochemical datasets, and is able to take the
biosynthetic/structural properties of the chemical compounds into
account while doing so.

## Installation

The developmental version of the package can be installed directly from
GitHub using the `install_github` function from the devtools package.

``` r
# Install devtools if not already installed
install.packages("devtools")

library("devtools")

install_github("hpetren/chemdiv")
```

## Usage

Detailed usage notes can be found with `help(chemdiv)` and in the
documentation for each function. In short, two types of data is
required. First, a dataset on the relative proportions of different
chemical compounds in different samples. Second, a dataset with the
name, SMILES and InChIKey of all the compounds found in the samples.

#### Compound classification and dissimilarities

Functions `NPCTable` and `compDis` classify chemical compounds and
calculate dissimilarities between them.

#### Diversity calculations

Functions `calcDiv`, `calcBetaDiv` and `calcDivProf` are used to
calculate different kinds of (alpha)-diversity and evenness (`calcDiv`),
beta-diversity (`calcBetaDiv`) and diversity profiles (`calcDivProf`)
for the samples.

#### Sample dissimilarities

Function `sampleDis` is used to calculate Generalized UniFrac
dissimilarities or Bray-Curtis dissimilarities between samples.

#### Molecular network and properties

Function `molNet` creates a molecular network of the chemical compounds
and calculates some network properties.

#### Chemodiversity and network plots

`molNetPlot` and `chemDivPlot` are used to conveniently create basic
plots of the molecular network and different types of chemodiversity
measurements created by various functions in the package.
