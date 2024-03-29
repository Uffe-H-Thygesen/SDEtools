% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdetools.R
\name{QuadraticVariation}
\alias{QuadraticVariation}
\title{Discretized quadratic variation of a stochastic process}
\usage{
QuadraticVariation(X)
}
\arguments{
\item{X}{numeric vector containing the process}
}
\value{
A numeric vector, same length as X, giving the discretized quadratic variation as a function of time
}
\description{
Discretized quadratic variation of a stochastic process
}
\examples{
## Quadratic Variation of Brownian motion
times <- seq(0,10,0.01)
B <- rBM(times)
plot(times,QuadraticVariation(B))

## Quadratic Variation of an Ito integral
G <- cos(B)
X <- itointegral(G,B)
plot(times,QuadraticVariation(X))
lines(times,itointegral(G^2,times))

}
