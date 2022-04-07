#' @title Prediction
#' @description \code{predGRR} This function provides predictive values for future observation
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr inset
#' @importFrom magrittr raise_to_power
#' @importFrom magrittr set_colnames
#' @importFrom purrr invoke
#' @importFrom purrr map
#' @param res the output of `GRR.MSC`, `pGRR.MSC`, `AM.pGRR`
#' @param newX a matrix of penalized explanatory variables for future observation
#' @param newZ a matrix of non-penalized explanatory variables, with intercept if need,
#'     for future observation
#' @return predictive values
#' @export
#' @examples
#' #predGRR(res, newX)

predGRR <- function(res, newX, newZ){

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
        OLS = res$mu + newX %*% res$coefOLS %>% drop,
        GRR = res$mu + newX %*% res$coefGRR %>% drop
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
    B1 <- map(1:k, ~cpoly(newX[,.x])) %>% invoke(cbind, .) %>% cbind(1, .) %>%
      set_colnames(NULL)

    TAU <- res$knots
    B2 <- map(1:k, ~ctrunc(newX[,.x], TAU[,.x])) %>% invoke(cbind, .) %>%
      set_colnames(NULL)

    return(drop(
      B1%*%res$gamma + B2%*%res$beta
    ))
  }
}
