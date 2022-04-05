#' @title Generalized ridge regression with the optimal ridge parameters based on
#'   model selection criterion minimization method
#' @description \code{GRR.MSC} This function provides generalized ridge estimator
#'   with the optimal ridge parameters based on model selection criterion minimization method
#'
#' @importFrom MASS ginv
#' @importFrom magrittr "%>%"
#' @importFrom magrittr extract
#' @importFrom magrittr inset
#' @importFrom magrittr set_names
#' @importFrom purrr map2_dbl
#' @param y a vector of a response variable
#' @param X a matrix of explanatory variables without intercept
#' @param MSC a model selection criterion; GCp or EGCV
#' @param alpha a value (>=2) expressing penalty strength for MSC
#' @param n sample size
#' @param tol tolerance for rank deficient
#' @return estimation results
#' @export
#' @examples
#' #GRR.MSC(y, X)

GRR.MSC <- function(y, X, MSC=c("GCp", "EGCV"), alpha=log(n), n=length(y), tol=1e-12){

  MSC <- MSC[1]
  cand <- NULL

  X0 <- X
  X <- scale(X0, scale=F)
  k <- ncol(X)

  X. <- t(X)
  M <- X. %*% X
  Minv <- ginv(M)
  X.y <- X. %*% y

  Xsvd <- svd(X)
  dsq <- Xsvd$d %>% extract(. > tol)
  d <- dsq^2
  m <- length(dsq)

  P1 <- Xsvd$u[,1:m]
  Q1 <- Xsvd$v[,1:m]

  mu <- mean(y)
  BetaLSE <- Minv %*% X.y %>% drop
  yLSE <- mu + X%*%BetaLSE %>% drop

  b <- (n - m - 1)/n

  if(b == 0)
  {
    sig0 <- 0
  } else
  {
    sig0 <- sum((y-yLSE)^2)/n
  }

  z <- t(P1) %*% y %>% drop
  z2 <- z^2
  t <- sort(z2)
  t0 <- c(0, t)

  c1 <- cumsum(t0)
  c2 <- sum(1/t) - c(0, cumsum(1/t))

  if(MSC == "GCp")
  {
    if(sig0 == 0){stop("GCp cannot be defined")}

    h <- alpha*sig0/(2*b)
  } #end if GCp

  if(MSC == "EGCV")
  {
    aa <- 0:(m-1)
    Rl <- t0[-(m+1)]; Rr <- t

    if(alpha == 2)
    {
      if(sig0 != 0)
      {
        sa <- (n*sig0 + c1[-(m+1)]) / (n - m - 1 + aa)
        a.ast <- which(Rl < sa & sa <= Rr)

        h <- sa[a.ast]
      } else if(b == 0)
      {
        h <- t[1]
      } else
      {
        h <- 0
      }
    } else
    {
      if(sig0 == 0 & b != 0)
      {
        h <- 0
      } else
      {
        siginf <- sum((y - mu)^2)/n

        anb1 <- aa + n*b
        anb2 <- (
          anb1^2 - alpha*(alpha-2)*c2[-(m+1)]*(n*sig0 + c1[-(m+1)])
        ) %>% inset(.<0, Inf) %>% sqrt
        Xi <- ( anb1 - anb2 ) / ( (alpha-2)*c2[-(m+1)] )

        A <- which(Rl < Xi & Xi <= Rr)
        S <- Xi[A]

        if(n*alpha*siginf > 2*(n-1)*t[m])
        {
          A <- c(A, m)
          S <- c(S, t[m])
        }

        if(length(S) == 1)
        {
          h <- S
        } else
        {
          phia <- function(h, a){
            c1a <- c1[a+1]
            c2a <- c2[a+1]

            out1 <- sig0 + (c1a + c2a*(h^2))/n
            out2 <- (b + (a + c2a*h)/n)^alpha

            return(out1 / out2)
          }

          h <- map2_dbl(S, A-1, phia) %>% which.min %>% S[.]
        }

        cand <- data.frame(
          value = S, a = A-1
        )
      }
    }
  } #end if EGCV

  v <- (1 - (h/z2)) %>% inset(.<0, 0)
  BetaGRR <- Q1 %*% diag(v) %*% t(Q1) %*% BetaLSE %>% drop
  yGRR <- mu + X%*%BetaGRR %>% drop

  R2 <- (
    1 - ( c(sum((y - yLSE)^2), sum((y - yGRR)^2)) / sum((y - mu)^2) )
  ) %>% set_names(c("LSE", "GRR"))

  return(
    list(
      mu = mu,
      coefLSE = BetaLSE,
      coefGRR = BetaGRR,
      predLSE = yLSE,
      predGRR = yGRR,
      theta = (
        d*h / (z2 - h)
      ) %>% inset(.<0, Inf) %>% c(numeric(k-m)),
      R2 = R2,
      svd = Xsvd,
      cand = cand
    )
  )
}
