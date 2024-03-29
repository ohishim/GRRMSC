#' @title Prediction (v0.2.0)
#' @description \code{predGRR} This function provides predictive values for future observation
#'
#' @importFrom magrittr "%>%" inset raise_to_power set_colnames
#' @importFrom purrr exec map
#' @param res the output of `GRR.MSC`, `pGRR.MSC`, `AM.pGRR`
#' @param newX a matrix of penalized explanatory variables for future observation
#' @param newZ a matrix of non-penalized explanatory variables, with intercept if need,
#'     for future observation
#' @return predictive values;
#'   if `class(res) == "GRR.MSC"`, this is a list object which has two predictive values
#'   obtained by OLS and GRR, respectively
#' @export
#' @examples
#' #predGRR(res, newX)

predGRR <- function(res, newX, newZ=NULL){

  if(class(res) == "GRR.MSC")
  {
    if(!is.null(res$Xcenter))
    {
      nn <- nrow(newX); k <- ncol(newX)
      newX0 <- newX
      newX <- newX0 - matrix(res$Xcenter, nn, k, byrow=T)
    }

    return(
      list(
        OLS = res$intercept + newX %*% res$coefOLS %>% drop,
        GRR = res$intercept + newX %*% res$coefGRR %>% drop
      )
    )
  }

  if(class(res) == "pGRR.MSC")
  {
    return(drop(
      newX %*% res$beta + newZ %*% res$gamma
    ))
  }

  if(class(res) == "AM.pGRR")
  {
    B1 <- map(1:k, ~cpoly(newX[,.x])) %>% exec(cbind, !!!.) %>% cbind(1, .) %>%
      set_colnames(NULL)

    TAU <- res$knots
    B2 <- map(1:k, ~ctrunc(newX[,.x], TAU[,.x])) %>% exec(cbind, !!!.) %>%
      set_colnames(NULL)

    return(drop(
      B1%*%res$gamma + B2%*%res$beta
    ))
  }
}
