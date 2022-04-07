
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
by a method proposed by Ohishi, Yanagihara & Wakaki (2020).

To predict for future observation, you can use `GRR.pred` like this:

``` r
# newX: a matrix of explanatory variables for future observation

pred <- GRR.pred(newX, res)
```

This package also has the following functions:

-   `pGRR.MSC`: for partial generalized ridge regression.  
-   `AM.pGRR`: for additive model via partial generalized ridge
    regression (Fukui <i>et al</i>., 2020).  
-   `FM`: for generation of simulation data via fixed model (Ohishi,
    Yanagihara & Fujikoshi, 2020).

## References

1.  Fukui, K., Ohishi, M., Yamamura, M. & Yanagihara, H. (2020). A fast
    optimization method for additive model via partial generalized ridge
    regression. <i>Smart Innov. Syst. Tech.</i>, <b>193</b>, 279-290.
    \[[link](https://doi.org/10.1007/978-981-15-5925-9_24)\].  
2.  Ohishi, M., Yanagihara, H. & Fujikoshi, Y. (2020). A fast algorithm
    for optimizing ridge parameters in a generalized ridge regression by
    minimizing a model selection criterion. <i>J. Statist. Plann.
    Inference</i>, <b>204</b>, 187-205.
    \[[link](https://doi.org/10.1016/j.jspi.2019.04.010)\]  
3.  Ohishi, M., Yanagihara, H. & Wakaki, H. (2020). Optimization of
    generalized *C*<sub>*p*</sub> criterion for selecting ridge
    parameters in generalized ridge regression. <i>Smart Innov. Syst.
    Tech.</i>, <b>193</b>, 267-278.
    \[[link](https://doi.org/10.1007/978-981-15-5925-9_23)\]
