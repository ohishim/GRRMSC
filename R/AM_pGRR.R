#' @title Additive model via partial generalized ridge regression
#' @description \code{AM.pGRR} This function estimates additive model via partial generalized ridge regression
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr inset
#' @importFrom magrittr raise_to_power
#' @importFrom magrittr set_colnames
#' @importFrom magrittr set_class
#' @importFrom purrr invoke
#' @importFrom purrr map
#' @importFrom stats quantile
#' @param y a vector of a response variable
#' @param X a matrix of explanatory variables without intercept
#' @param k the number of explanatory variables
#' @param q a vector or an integer which is the number of knots;
#'     if vector, the number of knots is different for each variable
#' @return a list object with "AM.pGRR" class which has the following elements:
#' \tabular{ll}{
#'   `B1`, `B2` \tab
#'     matrices of explanatory variables transformed by cubic truncated power basis function;
#'     `B2` is penalized and `B1` is not \cr \tab \cr
#'   `knots` \tab
#'     a matrix of knots \cr \tab \cr
#'   `beta` \tab
#'     estimates for `B2` \cr \tab \cr
#'   `gamma` \tab
#'     estimates for `B1` \cr \tab \cr
#'   `fitted.values` \tab
#'     fitted values \cr \tab \cr
#'   `R2` \tab
#'     the coefficient of determination \cr
#' }
#' @export
#' @examples
#' #AM.pGRR(y, X)

AM.pGRR <- function(y, X, k=ncol(X), q=k){

  if(length(q) == 1){q <- rep(q, k)}

  B1 <- map(1:k, ~cpoly(X[,.x])) %>% invoke(cbind, .) %>% cbind(1, .) %>%
    set_colnames(NULL)

  TAU <- map(1:k, ~knots(X[,.x], q[.x])) %>% invoke(cbind, .)
  B2 <- map(1:k, ~ctrunc(X[,.x], TAU[,.x])) %>% invoke(cbind, .) %>%
    set_colnames(NULL)

  out <- c(
    list(
      B1 = B1,
      B2 = B2,
      knots = TAU
    ),
    pGRR.MSC(y, B2, B1)
  ) %>% set_class("AM.pGRR")

  return(out)
}
