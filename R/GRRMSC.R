#' @title Generalized ridge regression optimized based on model selection criterion minimization method (v1.0.0)
#' @description \code{GRR.MSC} This function provides generalized ridge estimator
#'   with the optimal ridge parameters based on model selection criterion minimization method
#'
#' @importFrom MASS ginv
#' @importFrom magrittr "%>%" extract inset set_names set_class multiply_by
#' @importFrom purrr exec map map2_dbl
#' @param y a vector of a response variable
#' @param X a matrix of explanatory variables without intercept
#' @param MSC a model selection criterion:
#'     `"EGCV"` (default), `"GCV"`, `"GCp"`, `"Cp"`,`"MCp"`,
#'     `"GIC"`, `"AIC"`, `"HQC"`, `"BIC"`, `"AICc"` or `"GCVc"`
#' @param alpha a value (>=2) expressing penalty strength for `"EGCV"`, `"GCp"`,
#'     `"GIC"` or `"GCVc"`; default for `"EGCV"` and `"GIC"` is `log(n)`;
#'     default for `"GCp"` is `"optimized"`; default for `"GCVc"` is 1
#' @param n sample size
#' @param intercept if `FALSE`, intercept is removed from the model
#' @param centering if `X` is already centralized, set `FALSE`
#' @param tol tolerance for rank deficient
#' @return a list object with "GRR.MSC" class which has the following elements:
#' \item{intercept}{the estimate for intercept;
#'     if `centering=TRUE`, this is value for centralized explanatory variables}
#'
#' \item{coefOLS}{the ordinary least squares (OLS) estimates}
#'
#' \item{coefGRR}{the generalized ridge regression (GRR) estimates}
#'
#' \item{fitOLS}{the fitted values by OLS}
#'
#' \item{fitGRR}{the fitted values by GRR}
#'
#' \item{theta}{the optimal ridge parameters}
#'
#' \item{R2}{the coefficient of determination}
#'
#' \item{svd}{the output of `svd(X)`, where `X` is centralized if `centering=TRUE`}
#'
#' \item{cand}{candidates of `h` which gives the optimal ridge parameters}
#'
#' \item{MSC}{model selection criterion used to optimize ridge parameters}
#'
#' \item{alpha}{the value of penalty strength used in `MSC`}
#'
#' \item{Xcenter}{if `centering=TRUE`, a vector of means for each explanatory variable; if not, `NULL`}
#' @export
#' @examples
#' #GRR.MSC(y, X)

GRR.MSC <- function(
  y, X, MSC="EGCV", alpha="default", n=length(y),
  intercept=TRUE, centering=TRUE, tol=1e-12
){
  ##############################################################################
  ###   preparation
  ##############################################################################

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
  if(MSC == "AIC"){MSC <- "GIC"; alpha <- 2}
  if(MSC == "HQC"){MSC <- "GIC"; alpha <- 2*log(log(n))}
  if(MSC == "BIC"){MSC <- "GIC"; alpha <- log(n)}

  if(alpha == "default")
  {
    if(MSC %in% c("EGCV", "GIC")){alpha <- log(n)}
    if(MSC == "GCp"){alpha <- "optimized"}
    if(MSC == "AICc"){alpha <- "none"}
    if(MSC == "GCVc"){alpha <- 1}
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
  } else
  {#---   except GCp   ---------------------------------------------------------

    t <- sort(z2)
    t0 <- c(0, t)

    c1 <- cumsum(t0)
    c2 <- sum(1/t) - c(0, cumsum(1/t))

    aa <- 0:(m-1)
    Rl <- t0[-(m+1)]; Rr <- t

    if(MSC == "EGCV")
    {#---   EGCV   -------------------------------------------------------------

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

    if(MSC == "GIC")
    {#---   GIC   --------------------------------------------------------------

      if(sig0 == 0)
      {
        h <- 0
      } else
      {
        anb1 <- n
        anb2 <- (
          anb1^2 - (alpha^2)*c2[-(m+1)]*(n*sig0 + c1[-(m+1)])
        ) %>% inset(.<0, Inf) %>% sqrt
        Xi <- ( anb1 - anb2 ) / ( alpha*c2[-(m+1)] )

        A <- which(Rl < Xi & Xi <= Rr)
        S <- Xi[A]

        if(alpha*siginf > 2*t[m])
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
            out2 <- alpha*(1 + m - a - c2a*h)/n

            return(out1*exp(out2))
          }

          h <- map2_dbl(S, A-1, phia) %>% which.min %>% S[.]
        }

        cand <- data.frame(
          value = S, a = A-1
        )
      }
    } #end if GIC

    if(MSC == "AICc")
    {#---   AICc   -------------------------------------------------------------

      h0s <- (3 + m - aa - n)/c2[-(m+1)]
      a0 <- which(Rl < h0s & h0s <= Rr)

      if(length(a0) == 0)
      {
        exists.a0 <- FALSE
        a0 <- h0 <- 0
      } else
      {
        exists.a0 <- TRUE
        h0 <- h0s[a0]
        a0 <- a0 - 1
      }

      if(sig0 == 0 & exists.a0 == FALSE)
      {
        h <- 0
      } else
      {
        aa <- a0:(m-1); m0 <- length(aa)
        tt <- c(h0, t[aa+1])
        Rl <- tt[-(m0+1)]; Rr <- tt[-1]

        psi.coef <- cbind(
          c2[aa+1]^2,
          (n - 5 - 2*m + 2*aa)*c2[aa+1],
          (n - 3 - m + aa)^2,
          (1 - n)*(n*sig0 + c1[aa+1])
        )

        cand0 <- map(1:m0, ~{
          res <- solve3eq(psi.coef[.x,])
          if(res$discriminant == "positive")
          {
            return(cbind(range(res$real.root), aa[.x]))
          } else
          {
            return(c(res$real.root[1], aa[.x]))
          }
        }) %>% exec(rbind, !!!.) %>% as.data.frame %>% set_names(c("value", "a"))

        A <- which(Rl[cand0$a+1] < cand0$value & cand0$value <= Rr[cand0$a+1])
        S <- cand0$value[A]

        if(n*(n-1)*siginf > t[m]*((n-3)^2))
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
            out2 <- (2*(2 + m - a - c2a*h)) / (n - 3 - m + a + c2a*h)

            return(out1*exp(out2))
          }

          h <- map2_dbl(S, A-1, phia) %>% which.min %>% S[.]
        }

        cand <- data.frame(
          value = S, a = A-1
        )
      }
    } #end if AICc

    if(MSC == "GCVc")
    {#---   GCVc   -------------------------------------------------------------

      h0s <- (1 + m - aa + alpha - n)/c2[-(m+1)]
      a0 <- which(Rl < h0s & h0s <= Rr)

      if(length(a0) == 0)
      {
        exists.a0 <- FALSE
        a0 <- h0 <- 0
      } else
      {
        exists.a0 <- TRUE
        h0 <- h0s[a0]
        a0 <- a0 - 1
      }

      if(sig0 == 0 & exists.a0 == FALSE)
      {
        h <- 0
      } else if(n*siginf > (n - 1 - alpha)*t[m])
      {
        h <- t[m]
      } else
      {
        aa <- a0:(m-1)
        tt <- c(h0, t[aa+1])
        Rl <- tt[-length(tt)]; Rr <- tt[-1]

        sa <- (n*sig0 + c1[aa+1]) / (n - m - 1 + aa - alpha)
        a.ast <- which(Rl < sa & sa <= Rr)

        h <- sa[a.ast]
      }
    } #end if GCVc
  }

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
    intercept = mu,
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
