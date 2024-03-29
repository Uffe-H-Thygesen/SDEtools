% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdetools.R
\name{rvBM}
\alias{rvBM}
\title{Simulate a multivariate Brownian motion}
\usage{
rvBM(times, n = 1, sigma = rep(1, n), B0 = rep(0, n), u = rep(0, n))
}
\arguments{
\item{times}{numeric vector of time points where the Brownian motion should be simulated}

\item{n}{numeric scalar, dimension of Brownian motion}

\item{sigma}{noise level in the Brownian motion. A vector of length n, one for each dimension, or a scalar, in which case the same level is applied to all dimensions.}

\item{B0}{initial condition, applicable at time t=0. A vector of length n, one for each dimension, or a scalar, in which case the same initial condition is applied to all dimensions.}

\item{u}{Drift. An optional numeric n-vector which defaults to a vector of zeros. When supplied, a linear drift (bias) is added to each component.}
}
\value{
a numeric array, n*length(times), each column containing a sample path
}
\description{
Simulate a multivariate Brownian motion
}
\examples{
times <- 0:10
VB <- rvBM(times,3)
matplot(times,VB,type="b")

}
