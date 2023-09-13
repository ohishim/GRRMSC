#' @title Partial generalized ridge regression optimized based on model selection criterion minimization method (v0.2.0)
#' @description \code{pGRR.MSC} This function provides partial generalized ridge estimator
#'   with the optimal ridge parameters based on model selection criterion minimization method
#'
#' @importFrom magrittr "%>%" set_class
#' @param y a vector of a response variable
#' @param X a matrix of penalized explanatory variables without intercept
#' @param Z a matrix of non-penalized explanatory variables with intercept if need
#' @param MSC a model selection criterion; see the function `GRRMSC`
#' @param alpha a value (>=2) expressing penalty strength for `MSC`; see the function `GRRMSC`
#' @param n sample size
#' @param tol tolerance for rank deficient
#' @return a list object with "pGRR.MSC" class which has the following elements:
#' \item{beta}{estimates for `X`}
#'
#' \item{gamma}{estimates for `Z`}
#'
#' \item{fitted.values}{fitted values}
#'
#' \item{R2}{the coefficient of determination}
#' @export
#' @examples
#' #pGRR.MSC(y, X)

pGRR.MSC <- function(
  y, X, Z, MSC="EGCV", alpha="default", n=length(y), tol=1e-12
){
  Z. <- t(Z)
  W <- Z. %*% Z
  Winv <- solve(W)
  Z.y <- Z. %*% y
  Z.X <- Z. %*% X
  ZWinv <- Z %*% Winv

  res <- GRR.MSC(
    y - ZWinv%*%Z.y %>% drop,
    X - ZWinv%*%Z.X,
    MSC=MSC, alpha=alpha, n=n,
    intercept=FALSE, centering=FALSE, tol=tol
  )

  Beta <- res$coefGRR
  Gamma <- Winv %*% (Z.y - Z.X%*%Beta) %>% drop

  yh <- (X%*%Beta + Z%*%Gamma) %>% drop

  out <- list(
    beta = Beta,
    gamma = Gamma,
    fitted.values = yh,
    R2 = 1 - sum((y-yh)^2)/sum((y-mean(y))^2)
  ) %>% set_class("pGRR.MSC")

  return(out)
}
