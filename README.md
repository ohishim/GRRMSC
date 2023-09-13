
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GRRMSC (v1.2.0)

<!-- badges: start -->
<!-- badges: end -->

This package provides generalized ridge estimator with the optimal ridge
parameters based on a model selection criterion minimization method.

**cite this package**:  
Ohishi, M. (2023). GRRMSC: Generalized ridge regression optimized based
on model selection criterion minimization method. R package version
1.2.0. <https://github.com/ohishim/GRRMSC>

## Installation

You can install the R package GRRMSC from [GitHub](https://github.com/)
with:

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
with $\alpha=\log n$, where $n$ (and `n`) is sample size. To use other
criteria, you can set a parameter `MSC` at “GCV”, “GCp”, “Cp”, “MCp”,
“GIC”, “AIC”, “HQC”, “BIC”, “AICc” or “GCVc”. For `MSC="EGCV"`,
`MSC="GCp"`, `MSC="GIC"`, and `MSC="GCVc"`, `GRR.MSC` has a parameter
`alpha` to change a value of $\alpha$. When you set `MSC="GCp"` and
`alpha="optimized"`, a value of $\alpha$ in generalized $C_p$ criterion
is optimized by a method proposed by Ohishi, Yanagihara & Wakaki (2020).

To predict for future observation, you can use `predGRR` like this:

``` r
# newX: a matrix of explanatory variables for future observation

pred <- predGRR(res, newX)
```

This package also has the following functions:

- `pGRR.MSC`: for partial generalized ridge regression.  
- `AM.pGRR`: for additive model via partial generalized ridge regression
  (Fukui <i>et al</i>., 2020).  
- `FM`: for generation of simulation data via fixed model (Ohishi,
  Yanagihara & Fujikoshi, 2020).

Note that `predGRR` also performs for `pGRR.MSC` and `AM.pGRR`.

## References

1.  Fukui, K., Ohishi, M., Yamamura, M. & Yanagihara, H. (2020). A fast
    optimization method for additive model via partial generalized ridge
    regression. <i>Smart Innov. Syst. Tech.</i>, <b>193</b>, 279-290.
    doi:
    [10.1007/978-981-15-5925-9_24](https://doi.org/10.1007/978-981-15-5925-9_24)  
2.  Ohishi, M., Yanagihara, H. & Fujikoshi, Y. (2020). A fast algorithm
    for optimizing ridge parameters in a generalized ridge regression by
    minimizing a model selection criterion. <i>J. Statist. Plann.
    Inference</i>, <b>204</b>, 187-205. doi:
    [10.1016/j.jspi.2019.04.010](https://doi.org/10.1016/j.jspi.2019.04.010)  
3.  Ohishi, M., Yanagihara, H. & Wakaki, H. (2020). Optimization of
    generalized $C_p$ criterion for selecting ridge parameters in
    generalized ridge regression. <i>Smart Innov. Syst. Tech.</i>,
    <b>193</b>, 267-278. doi:
    [10.1007/978-981-15-5925-9_23](https://doi.org/10.1007/978-981-15-5925-9_23)
