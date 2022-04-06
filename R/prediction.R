#' @title Prediction
#' @description \code{GRR.pred} This function provides predictive values for future observation
#'
#' @importFrom magrittr "%>%"
#' @param newX a matrix of explanatory variables for future observation
#' @param res the output of GRR.MSC
#' @return predictive values
#' @export
#' @examples
#' #GRR.pred(newX, res)

GRR.pred <- function(newX, res){

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
