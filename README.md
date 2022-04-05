
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GRRMSC

<!-- badges: start -->
<!-- badges: end -->

This package provides generalized ridge estimator with the optimal ridge
parameters based on a model selection criterion minimization method.

## Installation

You can install the development version of GRRMSC from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ohishim/GRRMSC")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(GRRMSC)

# y: a vector of a response variable
# X: a matrix of explanatory variables without intercept

res <- GRR.MSC(y, X)
```

In default, ridge parameters are optimized by minimizing EGCV criterion
with *α* = log *n*, where *n* is sample size. To change the value of
*α*, `GRR.MSC` has a parameter `alpha`. You can use *G**C*<sub>*p*</sub>
criterion by setting a parameter of `MSC="GCp"`.
