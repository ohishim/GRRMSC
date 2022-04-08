#' @title Generation of simulation data
#' @description \code{FM} This function generates simulation data via fixed model
#'
#' @importFrom magrittr "%>%"
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @param n sample size
#' @param k the number of explanatory variables
#' @param rho correlation between explanatory variables
#' @param sig `sig^2` is an error variance
#' @param beta true regression coefficients (1 or 2)
#' @param k0 option when `beta=1`; the number of non-zero coefficients
#' @return a list object of simulation data which has the following elements:
#' \item{y}{a vector of a response variable}
#'
#' \item{X}{a matrix of explanatory variables}
#'
#' \item{beta}{a vector of the true coefficients}
#' @export
#' @examples
#' #FM(n, k)

FM <- function(n, k, rho=0.99, sig=1, beta=1, k0=floor(k/2)){

  PHI <- function(rho, k){
    Out <- diag(1, k, k)

    for(i in 1:(k-1))
    {
      for(j in (i+1):k)
      {
        Out[i, j] <- Out[j, i] <- rho^(j-i)
      }
    }

    return(Out)
  }

  matsq <- function(M){
    eigM <- eigen(M)
    return(
      eigM$vectors %*% (sqrt(eigM$values)*t(eigM$vectors))
    )
  }

  PSI <- PHI(rho, k)
  PSIsq <- matsq(PSI)

  X0 <- runif(n*k, -1, 1) %>% matrix(n, k)
  X <- X0 %*% PSIsq

  if(beta == 1)
  {
    Beta <- c(rep(1, k0), numeric(k-k0))
  } else if(beta == 2)
  {
    jj <- 1:k
    Beta <- sqrt( ( 12*n*(n-1) )/( 4*(n^2) + 6*n - 1 ) ) *
      ( (-1)^(jj-1) * (1 - (jj-1)/n - 1/(2*n)) )
  }

  Eta <- X %*% Beta %>% drop
  y <- Eta + rnorm(n, 0, sig)

  return(
    list(
      y = y,
      X = X,
      beta = Beta
    )
  )
}
