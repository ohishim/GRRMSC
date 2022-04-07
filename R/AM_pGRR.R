#' @title Additive model via partial generalized ridge regression
#' @description \code{AM.pGRR} This function estimates additive model via partial generalized ridge regression
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr inset
#' @importFrom magrittr raise_to_power
#' @importFrom purrr invoke
#' @importFrom purrr map
#' @param y a vector of a response variable
#' @param X a matrix of explanatory variables without intercept
#' @param k the number of explanatory variables
#' @param q a vector or an integer which is the number of knots;
#'     if vector, the number of knots is different for each variable
#' @return estimation results
#' @export
#' @examples
#' #AM.pGRR(y, X)

AM.pGRR <- function(y, X, k=ncol(X), q=k){

  if(length(q) == 1){q <- rep(q, k)}

  cpoly <- function(x){cbind(x, x^2, x^3)}

  ctrunc <- function(x, q){
    tau <- quantile(x, seq(0, 1, length=q+2)[-c(1, q+2)])
    return(
      map(tau, ~{
        (x - .x) %>% raise_to_power(3) %>% inset(.<0, 0)
      }) %>% invoke(cbind, .)
    )
  }

  B1 <- map(1:k, ~cpoly(X[,.x])) %>% invoke(cbind, .) %>% cbind(1, .)
  B2 <- map(1:k, ~ctrunc(X[,.x], q[.x])) %>% invoke(cbind, .)

  return(pGRR.MSC(y, B2, B1))
}
