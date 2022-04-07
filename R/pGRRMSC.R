#' @title Partial generalized ridge regression optimized based on model selection criterion minimization method
#' @description \code{pGRR.MSC} This function provides partial generalized ridge estimator
#'   with the optimal ridge parameters based on model selection criterion minimization method
#'
#' @importFrom magrittr "%>%"
#' @param y a vector of a response variable
#' @param X a matrix of penalized explanatory variables without intercept
#' @param Z a matrix of non-penalized explanatory variables with intercept if need
#' @param MSC a model selection criterion; `"EGCV"`, `"GCp"`, `"GCV"`, `"Cp"`, or `"MCp"`
#' @param alpha a value (>=2) expressing penalty strength for `MSC` (only when `MSC` is `"EGCV"` or `"GCp"`)
#' @param n sample size
#' @param tol tolerance for rank deficient
#' @return estimation results
#' @export
#' @examples
#' #pGRR.MSC(y, X)

pGRR.MSC <- function(
  y, X, Z, MSC=c("EGCV", "GCp", "GCV", "Cp", "MCp"), alpha=log(n), n=length(y), tol=1e-12
){
  MSC <- MSC[1]

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

  return(
    list(
      beta = Beta,
      gamma = Gamma,
      fitted.values = yh,
      r2 = 1 - sum((y-yh)^2)/sum((y-mean(y))^2)
    )
  )
}
