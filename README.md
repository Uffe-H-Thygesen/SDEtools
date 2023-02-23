
# SDEtools

<!-- badges: start -->
<!-- badges: end -->

SDEtools provides functions for analyzing models based on stochastic differential equations (SDE's). These analysis include simulation of sample paths, solving the Kolmogorov equations that govern transition probabilities, state estimation, and optimal control.

## Installation

You can install the development version of SDEtools from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Uffe-H-Thygesen/SDEtools")
```

## Example

This is a basic example which shows you how to simulate a sample path, and compute the transition probability, for a double well model:

``` r
library(SDEtools)
## basic example code
f <- function(x) x-x^3
g <- function(x) 0.5

tv <- seq(0,6,0.01)
sim <- heun(f,g,tv,0)

plot(tv,sim$X,type="l",ylim=c(-1.5,1.5))

xi <- seq(-1.5,1.5,0.01)
xc <- cell.centers(xi)

D <- function(x) 0.5*g(x)^2
G <- fvade(f,D,xi,'r')

rho <- StationaryDistribution(G)

lines(200*rho,xc)
```

