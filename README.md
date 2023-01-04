
<!-- README.md is generated from README.Rmd. Please edit that file -->

# chemodiv

<!-- badges: start -->
<!-- badges: end -->

`chemodiv` is an R package for analysing chemodiversity of phytochemical
data

The package can be used to calculate various types of diversities and
dissimilarities for phytochemical datasets. This includes the use of
diversity indices and dissimilarity indices that incorporate the
biosynthetic and/or structural properties of the phytochemical compounds
for the calculations, resulting in comprehensive measures of
phytochemical diversity. A complete description of the package is
available in Petrén et al. 2023.

## Installation

The current version of the package can be installed from CRAN.
Alternatively, the developmental version of the package can be installed
from GitHub using the `install_github()` function from the `devtools`
package.

``` r
# Install current version
install.packages("chemodiv")

# Install developmental version
install.packages("devtools") # Install devtools if not already installed
library("devtools")
install_github("hpetren/chemodiv")
```

## Usage

Detailed usage notes can be found with `help(chemodiv)` and in the
documentation for each function. See the vignette for a worked example.
In short, two datasets are required. First, a dataset on the relative
relative abundance/concentration (i.e. proportion) of different
phytochemical compounds in different samples. Second, a dataset with the
compound name, SMILES and InChIKey for all the compounds in the first
dataset. The following functions can then be used:

#### Data formatting checks

Function `chemoDivCheck()` checks so that the datasets used by functions
in the package are correctly formatted.

#### Compound classification and dissimilarities

Function `NPCTable()` uses the deep-learning tool *NPClassifier* to
classify chemical compounds into groups largely corresponding to
biosynthetic pathways. The function `compDis()` calculates and outputs a
list of dissimilarity matrices with pairwise dissimilarities between
chemical compounds, based on their biosynthetic and/or structural
properties.

#### Diversity calculations

Functions `calcDiv()`, `calcBetaDiv()` and `calcDivProf()` are used to
calculate phytochemical diversity in different ways, using both
traditional indices and Hill numbers. `calcDiv()` calculates alpha
diversity and evenness. `calcBetaDiv()` calculates beta diversity.
`calcDivProf()` generates diversity profiles. Calculations of functional
Hill numbers utilize a dissimilarity matrix generated by the `compDis()`
function.

#### Sample dissimilarities

Function `sampDis()` is used to calculate Generalized UniFrac
dissimilarities or Bray-Curtis dissimilarities between samples.
Calculations of Generalized UniFrac dissimilarities utilizes
dissimilarity matrix generated by the `compDis()` function.

#### Molecular network and properties

Function `molNet()` creates a molecular network of the chemical
compounds and calculates some network properties using matrices
generated by the `compDis()` function.

#### Chemodiversity and network plots

`molNetPlot()` and `chemoDivPlot()` are used to conveniently create
basic plots of the molecular network and different types phytochemical
diversity and dissimilarity calculated by the other functions in the
package.

#### Shortcut function

Function `quickChemoDiv()` uses many of the other functions in the
package to in one simple step calculate and visualize chemodiversity for
users wanting to quickly explore their data using standard parameters.

## References

Petrén H., T.G. Köllner and R.R. Junker. 2023. Quantifying
chemodiversity considering biochemical and structural properties of
compounds with the R package *chemodiv*. New Phytologist doi:
10.1111/nph.186850.
