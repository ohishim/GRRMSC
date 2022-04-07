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
  }) %>% invoke(cbind, .)
}

usethis::use_data(cpoly, knots, ctrunc, internal=T, overwrite = T)

