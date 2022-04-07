#' @title Generalized ridge regression optimized based on model selection criterion minimization method
#' @description \code{GRR.MSC} This function provides generalized ridge estimator
#'   with the optimal ridge parameters based on model selection criterion minimization method
#'
#' @importFrom MASS ginv
#' @importFrom magrittr "%>%"
#' @importFrom magrittr extract
#' @importFrom magrittr inset
#' @importFrom magrittr set_names
#' @importFrom magrittr set_class
#' @importFrom magrittr multiply_by
#' @importFrom purrr map2_dbl
#' @param y a vector of a response variable
#' @param X a matrix of explanatory variables without intercept
#' @param MSC a model selection criterion; `"EGCV"`, `"GCp"`, `"GCV"`, `"Cp"`, or `"MCp"`
#' @param alpha a value (>=2) expressing penalty strength for `MSC` (only when `MSC` is `"EGCV"` or `"GCp"`)
#' @param n sample size
#' @param intercept if `FALSE`, intercept is removed from the model
#' @param centering if `TRUE`, `X` is centralized
#' @param tol tolerance for rank deficient
#' @return estimation results
#' @export
#' @examples
#' #GRR.MSC(y, X)

GRR.MSC <- function(
  y, X, MSC=c("EGCV", "GCp", "GCV", "Cp", "MCp"), alpha=log(n), n=length(y),
  intercept=TRUE, centering=TRUE, tol=1e-12
){
  ##############################################################################
  ###   preparation
  ##############################################################################

  MSC <- MSC[1]
  cand <- NULL

  if(centering)
  {
    X0 <- X
    X <- scale(X0, scale=F)
    Xcenter <- attributes(X)$`scaled:center`
  } else
  {
    Xcenter <- NULL
  }

  k <- ncol(X)

  if(MSC == "GCV"){MSC <- "EGCV"; alpha <- 2}
  if(MSC == "Cp"){MSC <- "GCp"; alpha <- 2}
  if(MSC == "MCp")
  {
    if(n-k-3 <= 0)
    {
      stop("MCp cannot be defined")
    } else
    {
      MSC <- "GCp"; alpha <- 2 + ( 4/(n-k-3) )
    }
  }

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

  ybar <- mean(y)
  mu <- ifelse(intercept, ybar, 0)
  BetaOLS <- Minv %*% X.y %>% drop
  yOLS <- mu + X%*%BetaOLS %>% drop

  b <- (n - m - 1)/n

  if(b == 0)
  {
    sig0 <- 0
  } else
  {
    sig0 <- sum((y-yOLS)^2)/n
  }

  siginf <- sum((y - ybar)^2)/n

  z <- t(P1) %*% y %>% drop
  z2 <- z^2

  ##############################################################################
  ###   ridge parameters optimization
  ##############################################################################

  if(MSC == "GCp")
  {#---   GCp   ----------------------------------------------------------------

    if(sig0 == 0){stop("GCp cannot be defined")}

    s2 <- sig0 / b

    if(alpha == "optimized")
    {
      if(k >= n-3){stop("alpha cannot be optimized")}

      t <- sort(z2/s2)
      t0 <- c(0, t)

      c1 <- cumsum(t0)
      c2 <- sum(1/t) - c(0, cumsum(1/t))

      q <- 1 - ( 2/(n-k-1) )

      Phia <- q*c2*(t0^2) + 2*c2*t0 - 2*(0:k) + q*c1
      alpha <- which.min(Phia) %>% t0[.] %>% multiply_by(2)
    }

    h <- alpha*s2/2
  } #end if GCp

  if(MSC == "EGCV")
  {#---   EGCV   ---------------------------------------------------------------

    t <- sort(z2)
    t0 <- c(0, t)

    c1 <- cumsum(t0)
    c2 <- sum(1/t) - c(0, cumsum(1/t))

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

  ##############################################################################
  ###   output
  ##############################################################################

  v <- (1 - (h/z2)) %>% inset(.<0, 0)
  BetaGRR <- Q1 %*% (v*t(Q1)) %*% BetaOLS %>% drop
  yGRR <- mu + X%*%BetaGRR %>% drop

  R2 <- (
    1 - ( c(sum((y - yOLS)^2), sum((y - yGRR)^2)) / (n*siginf) )
  ) %>% set_names(c("OLS", "GRR"))

  out <- list(
    mu = mu,
    coefOLS = BetaOLS,
    coefGRR = BetaGRR,
    fitOLS = yOLS,
    fitGRR = yGRR,
    theta = (
      d*h / (z2 - h)
    ) %>% inset(.<0, Inf) %>% c(numeric(k-m)),
    R2 = R2,
    svd = Xsvd,
    cand = cand,
    MSC = MSC,
    alpha = alpha,
    Xcenter = Xcenter
  ) %>% set_class("GRR.MSC")

  return(out)
}
