% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generators.R
\name{prob2pdf}
\alias{prob2pdf}
\title{Convert cell probabilities to (average) probability densities}
\usage{
prob2pdf(phi, xgrid, ygrid = c(0, 1))
}
\arguments{
\item{phi}{Vector of probabilities}

\item{xgrid}{Spatial grid}

\item{ygrid}{Optional grid in y-direction (see details)}
}
\description{
prob2pdf takes a vector of probabilities and a spatial grid (given by end-points and interfaces)
and computes the average p.d.f. in each grid cell
}
\details{
length(phi) must equal length(xgrid)-1, if ygrid is not given, or
(length(xgrid)-1)*(length(ygrid)-1) if ygrid is given

xgrid (and ygrid, if given) must be strictly increasing
}
\examples{
xgrid <- seq(-2,2,length=10)^3
phi <- diff(pnorm(xgrid))
f <- prob2pdf(phi,xgrid)
plot(utils::head(xgrid,-1),f,type="S")

}
