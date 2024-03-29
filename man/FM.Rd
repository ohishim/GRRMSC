% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FixedModel.R
\name{FM}
\alias{FM}
\title{Generation of simulation data (v0.1.0)}
\usage{
FM(n, k, rho = 0.99, sig = 1, beta = 1, k0 = floor(k/2))
}
\arguments{
\item{n}{sample size}

\item{k}{the number of explanatory variables}

\item{rho}{correlation between explanatory variables}

\item{sig}{\code{sig^2} is an error variance}

\item{beta}{true regression coefficients (1 or 2)}

\item{k0}{option when \code{beta=1}; the number of non-zero coefficients}
}
\value{
a list object of simulation data which has the following elements:
\item{y}{a vector of a response variable}

\item{X}{a matrix of explanatory variables}

\item{beta}{a vector of the true coefficients}
}
\description{
\code{FM} This function generates simulation data via fixed model
}
\examples{
#FM(n, k)
}
