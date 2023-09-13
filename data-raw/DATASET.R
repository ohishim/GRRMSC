## code to prepare `DATASET` dataset goes here

# usethis::use_data(DATASET, overwrite = TRUE)

#---   functions   -------------------------------------------------------------

cpoly <- function(x){cbind(x, x^2, x^3)}

knots <- function(x, q){
  quantile(x, seq(0, 1, length=q+2)[-c(1, q+2)])
}

ctrunc <- function(x, tau){
  map(tau, ~{
    (x - .x) %>% raise_to_power(3) %>% inset(.<0, 0)
  }) %>% exec(cbind, !!!.)
}

solve3eq <- function(X){

  a <- X[1]; b <- X[2]; c <- X[3]; d <- X[4]
  D <- -4*a*(c^3) - 27*(a^2)*(d^2) + (b^2)*(c^2) + 18*a*b*c*d - 4*(b^3)*d

  p <- (c/(3*a)) - ((b^2)/(9*(a^2)))
  q <- ((b^3)/(27*(a^3))) - (b*c/(6*(a^2))) + (d/(2*a))

  r0 <- (q^2) + (p^3) #+ 0*1i
  r00 <- ifelse(r0 >= 0, sqrt(r0), sqrt(-r0)*1i)
  r1 <- -q + r00
  r2 <- -q - r00

  omega1 <- (-1+sqrt(3)*1i)/2
  omega2 <- (-1-sqrt(3)*1i)/2

  if(r0 < 0)
  {
    s1 <- r1^(1/3)
    s2 <- r2^(1/3)
  } else
  {
    s1 <- ifelse(r1 > 0, r1^(1/3), -(-r1)^(1/3))
    s2 <- ifelse(r2 > 0, r2^(1/3), -(-r2)^(1/3))
  }

  h1 <- s1 + s2
  h2 <- omega1*s1 + omega2*s2
  h3 <- omega2*s1 + omega1*s2

  x1 <- h1 - (b/(3*a))
  x2 <- h2 - (b/(3*a))
  x3 <- h3 - (b/(3*a))

  if(D > 0)
  {
    type <- "positive"
    x <- Re(c(x1, x2, x3))
    x0 <- NULL
  } else if(D == 0)
  {
    type <- "zero"

    if(p == 0 & q == 0)
    {
      x <- rep(x1, 3)
    } else
    {
      x <- Re(c(x1, x2, x3))
    }

    x0 <- NULL
  } else
  {
    type <- "negative"
    x <- x1
    x0 <- c(x2, x3)
  }

  return(list(
    real.root = x,
    complex.root = x0,
    discriminant = type
  ))
}

usethis::use_data(cpoly, knots, ctrunc, solve3eq, internal=T, overwrite = T)

