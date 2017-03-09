library(SDEtools)

test_that("SDE simulation: Dimensions of states and noise channels",{
  times <- 0:10

  nx <- 3
  f <- function(t,x) -x
  g <- function(t,x) x
  sol <- euler(f,g,times,rep(1,nx))

  expect_equal(dim(sol$X),c(length(times),nx))

  nb <- 2
  g <- function(t,x) matrix(runif(nx*nb),ncol=nb,nrow=nx)
  sol <- euler(f,g,times,rep(1,nx))

  expect_equal(dim(sol$X),c(length(times),nx))
})

test_that("Linear SDEs",{
  sol <- dLinSDE(-1,1,1,S0=0.5)
  expect_equal(as.numeric(sol$eAt),exp(-1))
  expect_equal(as.numeric(sol$St),0.5)

  A <- matrix(c(0,-1,1,0),nrow=2,ncol=2)
  G <- diag(c(1,1))
  sol <- dLinSDE(A,G,2*pi,x0=c(1,0))
  expect_equal(as.numeric(sol$EX),c(1,0))
  expect_equal(sol$St,diag(rep(2*pi,2)))

  sol <- dLinSDE(-1,1,2,x0=1,u=1)
  expect_equal(as.numeric(sol$EX),1)

  sol <- dLinSDE(0,1,2,x0=1,u=matrix(c(1,3),nrow=1,ncol=2))
  expect_equal(as.numeric(sol$EX),5)
})
