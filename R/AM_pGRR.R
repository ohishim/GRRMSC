#' @title Additive model via partial generalized ridge regression (v0.3.0)
#' @description \code{AM.pGRR} This function estimates additive model via partial generalized ridge regression
#'
#' @importFrom magrittr "%>%" inset raise_to_power set_colnames set_class
#' @importFrom purrr exec map
#' @importFrom stats quantile
#' @param y a vector of a response variable
#' @param X a matrix of explanatory variables without intercept
#' @param k the number of explanatory variables
#' @param q a vector or an integer which is the number of knots;
#'     if vector, the number of knots is different for each variable
#' @param MSC a model selection criterion; see the function `GRRMSC`
#' @param alpha a value (>=2) expressing penalty strength for `MSC`; see the function `GRRMSC`
#' @param tol tolerance for rank deficient
#' @return a list object with "AM.pGRR" class which has the following elements:
#' \item{B1, B2}{matrices of explanatory variables transformed by cubic truncated power basis function;
#'     `B2` is penalized and `B1` is not}
#'
#' \item{knots}{a matrix of knots}
#'
#' \item{beta}{estimates for `B2`}
#'
#' \item{gamma}{estimates for `B1`}
#'
#' \item{fitted.values}{fitted values}
#'
#' \item{R2}{the coefficient of determination}
#' @export
#' @examples
#' #AM.pGRR(y, X)

AM.pGRR <- function(
    y, X, k=ncol(X), q=k,
    MSC="EGCV", alpha="default", tol=1e-12
){

  if(length(q) == 1){q <- rep(q, k)}

  B1 <- map(1:k, ~cpoly(X[,.x])) %>% exec(cbind, !!!.) %>% cbind(1, .) %>%
    set_colnames(NULL)

  TAU <- map(1:k, ~knots(X[,.x], q[.x])) %>% exec(cbind, !!!.)
  B2 <- map(1:k, ~ctrunc(X[,.x], TAU[,.x])) %>% exec(cbind, !!!.) %>%
    set_colnames(NULL)

  out <- c(
    list(
      B1 = B1,
      B2 = B2,
      knots = TAU
    ),
    pGRR.MSC(y, B2, B1, MSC=MSC, alpha=alpha, tol=tol)
  ) %>% set_class("AM.pGRR")
  return(out)
}
