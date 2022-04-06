
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
with *α* = log *n*, by setting parameters at `MSC="EGCV"` and
`alpha=log(n)` , where *n* (and `n`) is sample size. To use other
criteria, you can set a parameter `MSC` at “GCp”, “GCV”, “Cp”, or “MCp”.
For `MSC="EGCV"` or `MSC="GCp"`, `GRR.MSC` has a parameter `alpha` to
change a value of *α*. When you set `MSC="GCp"` and `alpha="optimized"`,
a value of *α* in generalized *C*<sub>*p*</sub> criterion is optimized
by a method proposed by Ohishi, Yanagihara & Wakaki (2020, KES-IDT-20;
[10.1007/978-981-15-5925-9\_23](https://doi.org/10.1007/978-981-15-5925-9_23)).
