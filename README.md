
<!-- README.md is generated from README.Rmd. Please edit that file -->

# chemdiv

<!-- badges: start -->
<!-- badges: end -->

chemdiv is an R package used to analyse phytochemical data. **The
package is in development and may not be fully functional.**

The package can be used to calculate various types of diversities and
dissimilarities for phytochemical datasets. This includes the use of
diversity indices and dissimilarity indices that incorporate the
biosynthetic and/or structural properties of the phytochemical compounds
for the calculations, resulting in comprehensive measures of
phytochemical diversity.

## Installation

The developmental version of the package can be installed directly from
GitHub in R using the `install_github` function from the `devtools`
package.

``` r
# Install devtools if not already installed
install.packages("devtools")

library("devtools")

install_github("hpetren/chemdiv")
```

## Usage

Detailed usage notes can be found with `help(chemdiv)` and in the
documentation for each function. In short, datasets are required. First,
a dataset on the proportions of phytochemical compounds in different
samples. Second, a dataset with the compound name, SMILES and InChIKey
for all the compounds in the first dataset.

#### Data formatting checks

Function `chemDivCheck` can be used to check so that the datasets used
by functions in the package are correctly formatted.

#### Compound classification and dissimilarities

Function `NPCTable` uses NPClassifier to classify chemical compounds
into groups largely corresponding to biosynthetic pathways. Function
`compDis` calculates and outputs a list of dissimilarity matrices with
pairwise dissimilarities between chemical compounds, based on their
biosynthetic and/or structural properties.

#### Diversity calculations

Functions `calcDiv`, `calcBetaDiv` and `calcDivProf` are used to
calculate phytochemical diversity in different ways, using both
traditional indices and Hill numbers. `calcDiv` calculates measures
alpha-diversity and evenness. `calcBetaDiv` calculates beta-diversity.
`calcDivProf` generates diversity profiles. Calculations of functional
Hill numbers utilize dissimilarity matrices generated by the `compDis`
function.

#### Sample dissimilarities

Function `sampleDis` is used to calculate Generalized UniFrac
dissimilarities or Bray-Curtis dissimilarities between samples.
Calculations of Generalized UniFrac dissimilarities utilize
dissimilarity matrices generated by the `compDis` function.

#### Molecular network and properties

Function `molNet` creates a molecular network of the chemical compounds
and calculates some network properties using matrices generated by the
`compDis` function.

#### Chemodiversity and network plots

`molNetPlot` and `chemDivPlot` are used to conveniently create basic
plots of the molecular network and different types phytochemical
diversity and dissimilarity calculated by the other functions in the
package.

#### Shortcut function

Function `quickChemDiv` uses many of the other functions in the package
to in one simple step visualize chemodiversity for users wanting to
quickly explore their data using standard parameters.

## Vignette

Under construction.
