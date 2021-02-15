
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BestTransform: Identify the transformation that best fits all the continous vectors in a given dataframe.

The `BestTransform` R package was designed to help find a transformation
that best fits all the continuous vectors in a given dataframe. There
are many techniques that have been developed in this aim, however each
has been subject to their own strengths/weaknesses, and it is unclear on
how to decide which will work best until the data is observed. This
package will look at a range of possible transformations and return the
best one, based on the user specified metric, wich can be chosen from
Minimum skewness, Shapiro P value and Pearson P value.

## Installation

You can install the most recent (devel) version of BestTransform from
github with:

``` r
# install.packages("devtools")
devtools::install_github("SrinivasaSaketh/BestTransform")
```

Or, you can download it from CRAN with:

``` r
install.packages("BestTransform")
```

## Example

``` r
library(BestTransform)
```

``` r
#Load the dataset
data <- iris
#Transform the data
trans <- BestTransform::BestTransform(data, "Species")  
#> Loading required package: lattice
#> Loading required package: ggplot2
```
