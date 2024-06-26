#' Simulate a sample path of Brownian motion
#'
#' @name rBM
#' 
#' @param times numeric vector of time points where the Brownian motion should be simulated
#' @param sigma numeric scalar, noise level in the Brownian motion, defaults to 1
#' @param B0 numeric scalar, initial condition (at time=0) for the Brownian motion
#' @param u numeric scalar, drift
#' @return numeric vector, of same length as times, containing the simulated value of the Brownian motion
#' @seealso \code{\link{rvBM}} for a vectorized version
#' @examples
#' times <- 0:10
#' plot(times,rBM(times),type="b",xlab="Time",ylab="Brownian motion")
#' 
#' @export
rBM <- function(times,sigma=1,B0=0,u=0)
{
  dt <- c(times[1],diff(times))
  dB <- rnorm(length(times),mean=u*dt,sd=sigma*sqrt(dt))
  B <- B0+cumsum(dB)
  return(B)
}

#' Simulate a multivariate Brownian motion
#'
#' @name rvBM
#' 
#' @param times numeric vector of time points where the Brownian motion should be simulated
#' @param n numeric scalar, dimension of Brownian motion
#' @param sigma noise level in the Brownian motion. A vector of length n, one for each dimension, or a scalar, in which case the same level is applied to all dimensions.
#' @param B0 initial condition, applicable at time t=0. A vector of length n, one for each dimension, or a scalar, in which case the same initial condition is applied to all dimensions.
#' @param u Drift. An optional numeric n-vector which defaults to a vector of zeros. When supplied, a linear drift (bias) is added to each component.
#' @return a numeric array, n*length(times), each column containing a sample path
#' @examples
#' times <- 0:10
#' VB <- rvBM(times,3)
#' matplot(times,VB,type="b")
#' 
#' @export
rvBM <- function(times,n=1,sigma=rep(1,n),B0=rep(0,n),u=rep(0,n))
{
  if(length(sigma)==1)sigma <- rep(sigma,n)
  if(length(B0)==1)B0 <- rep(B0,n)
  if(length(u)==1)u <- rep(u,n)
  return(sapply(1:n,function(i)rBM(times,sigma=sigma[i],B0=B0[i],u=u[i])))
}

#' Simulate a Brownian bridge
#'
#' @name rBrownianBridge
#' 
#' @param times numeric vector of time points where the Brownian motion should be simulated
#' @param sigma scalar; the noise level in the Brownian motion.
#' @param B0T initial and terminal condition; the fixed value of the Brownian motion at the first and last time points
#' @return a numeric vector, length(times), containing a sample path
#' @examples
#' times <- 10:20
#' BB <- rBrownianBridge(times,c(10,-10),2)
#' plot(times,BB,type="b")
#' 
#' @export
rBrownianBridge <- function(times,B0T=c(0,0),sigma=1)
{
    B <- rBM(times,sigma=sigma)
    n <- length(times)
    t1 <- times[1]
    t2 <- times[n]

    slope <- (B0T[2] - B0T[1] - B[n] + B[1])/(times[n]-times[1])

    BB <- B + B0T[1] - B[1] + slope * (times - times[1])
    return(BB)
}

#' Integrate one sample path of a stochastic process w.r.t. another,
#' returning the Ito integral, the Stratonovich integral, or the
#' "right hand rule".
#'
#' @name stochint
#' 
#' @param f numeric vector containing the integrator
#' @param g numeric vector, same length as f, containing the integrand
#' @param rule character vector indicating rule(s). Valid choices are "l", "r", or "c" for "left", "right", and "center", as well as combinations, e.g. c("l","c"). Optional and defaults to "l".
#' @return A numeric vector, same length as f, giving the "running integral", i.e. the integral as a function of the upper limit.
#' @examples
#' ## Integrating a cosine w.r.t. a sine
#' times <- seq(0,2*pi,length=21)
#' I <- stochint(cos(times),sin(times))
#' Ia <- 0.5*times+0.25*sin(2*times)  # Analytical result
#' matplot(times,I,type="l")
#' lines(times,Ia,col="blue",lwd=2)
#'
#' ## Integration of Brownian motion w.r.t. itself
#' times <- seq(0,10,0.01)
#' BM <- rBM(times)
#' I <- stochint(BM,BM,c("l","c","r"))
#' matplot(times,cbind(I$l,0.5*BM^2-0.5*times),type="l",xlab="Time",ylab="Left integral (Ito)",
#'          main="Integral of B.M. w.r.t itself")
#' matplot(times,cbind(I$r,0.5*BM^2+0.5*times),type="l",xlab="Time",ylab="Right integral",
#'          main="Integral of B.M. w.r.t itself")
#' matplot(times,cbind(I$c,0.5*BM^2),type="l",xlab="Time",ylab="Central integral (Stratonovich)",
#'          main="Integral of B.M. w.r.t itself")
#' 
#' @export
stochint <- function(f,g,rule="l")
{
  n <- length(f)
  dg <- diff(g)
  l <- c(0,cumsum(f[-n]*dg))
  r <- c(0,cumsum(f[-1]*dg))
  c <- 0.5*(l+r)
  return(data.frame(l=l,r=r,c=c)[,rule])
}

#' Discretized quadratic variation of a stochastic process 
#'
#' @name QuadraticVariation
#' 
#' @param X numeric vector containing the process
#' @return A numeric vector, same length as X, giving the discretized quadratic variation as a function of time
#' @examples
#' ## Quadratic Variation of Brownian motion
#' times <- seq(0,10,0.01)
#' B <- rBM(times)
#' plot(times,QuadraticVariation(B))
#'
#' ## Quadratic Variation of an Ito integral
#' G <- cos(B)
#' X <- itointegral(G,B)
#' plot(times,QuadraticVariation(X))
#' lines(times,itointegral(G^2,times))
#' 
#' @export
QuadraticVariation <- function(X)
{
  return(c(0,cumsum(diff(X)^2)))
}

#' Ito integral of a stochastic process w.r.t. another
#'
#' @name itointegral
#' 
#' @param G numeric vector containing the integrator
#' @param B numeric vector, same length as G, containing the integrand
#' @return A numeric vector, same length as G, giving the "running integral", i.e. the Ito integral as a function of the upper limit.
#' @examples
#' ## Integration of Brownian motion w.r.t. itself
#' times <- seq(0,10,0.01)
#' BM <- rBM(times)
#' I <- itointegral(BM,BM)
#' matplot(times,cbind(I,0.5*BM^2-0.5*times),type="l",xlab="Time",ylab="Ito integral",
#'          main="Integral of B.M. w.r.t itself")
#' 
#' @export
itointegral <- function(G,B)
{
  return(c(0,cumsum(G[-length(G)]*diff(B))))
}

#' Cross-variation of two stochastic processes
#'
#' @name CrossVariation
#'
#' @param X numeric vector containing the one process
#' @param Y numeric vector, same length as X, containing the other process
#' @return A numeric vector, same length as X, giving the quadratic cross-varation as function of time
#' @examples
#' ## Quadratic variation of Brownian motion
#' times <- seq(0,10,0.01)
#' B <- rBM(times)
#' plot(times,CrossVariation(B,B),type="l")
#' abline(0,1,lty="dashed")
#'
#' ## Verifying Ito's formula
#' X <- -times + B
#' h <- dhdx <- dh2dx2 <- function(x) exp(x)
#' Y <- h(X)
#' plot(times,Y,type="l")
#'
#' Yi <- h(X[1]) + itointegral(dhdx(X),X) + 0.5*itointegral(dh2dx2(X),CrossVariation(X,X))
#' lines(times,Yi,col="blue")
#' 
#' @export
CrossVariation <- function(X,Y)
{
  return(c(0,cumsum(diff(X)*diff(Y))))
}

#' Heun's method for simulation of a Stratonovich stochastic differential equation dX = f dt + g dB
#'
#' @name heun
#'
#' @param f function f(x) or f(t,x) giving the drift term in the SDE
#' @param g function g(x) or g(t,x) giving the noise term in the SDE
#' @param times numeric vector of time points, increasing
#' @param x0 numeric vector giving the initial condition
#' @param B sample path of (multivariate) Browniań motion. If omitted, a sample is drawn using rvBM
#' @param p Optional projection function that at each time poins projects the state, e.g. to enforce non-negativeness
#' @param h Optional stopping function; the simulation terminates if h(t,x)<=0
#' @examples
#' times <- seq(0,10,0.1)
#'
#' # Plot a sample path of the solution of the SDE dX=-X dt + dB
#' plot(times,heun(function(x)-x,function(x)1,times,0)$X,type="l")
#'
#' # Plot a sample path of the solution of the 2D time-varying SDE dX=(AX + FU) dt + G dB
#'
#' A <- array(c(0,-1,1,-0.1),c(2,2))
#' F <- G <- array(c(0,1),c(2,1))
#'
#' f <- function(t,x) A %*% x + F * sin(2*t)
#' g <- function(x) G
#'
#' BM <- rBM(times)
#'
#' plot(times,heun(f,g,times,c(0,0),BM)$X[,1],type="l")
#'
#' # Add solution corresponding to different initial condition, same noise
#'
#' lines(times,euler(f,g,times,c(1,0),BM)$X[,1],type="l",lty="dashed")
#'
#' # Solve geometric Brownian motion as a Stratonovich equation
#' times <- seq(0,10,0.01)
#' BM <- rBM(times)
#' plot(times,heun(function(x)-0.1*x,function(x)x,times,1,BM)$X,type="l",xlab="Time",ylab="X")
#'
#' # Add solution of corresponding Ito equation
#' lines(times,euler(function(x)0.4*x,function(x)x,times,1,BM)$X,lty="dashed")
#'
#' # Testing a case with non-commutative noise
#' BM2 <- rvBM(times,2)
#' matplot(times,heun(function(x)0*x,
#'   function(x)diag(c(1,x[1])),times,c(0,0),BM2)$X,xlab="Time",ylab="X",type="l")
#'
#' # Add solution to the corresponding Ito equation
#' matplot(times,euler(function(x)0*x,
#'   function(x)diag(c(1,x[1])),times,c(0,0),BM2)$X,xlab="Time",ylab="X",type="l",add=TRUE)
#' 
#' @export
heun <- function(f,g,times,x0,B=NULL,p=function(x)x,h=NULL,r=NULL,S0=1,Stau=runif(1))
{
  nx <- length(x0)
  nt <- length(times)

  ## Check if f is specified as a function(x) only, then convert to a function(t,x)
  if(length(formals(f))==1){
    ff <- function(t,x)f(x)
  }else
  {
    ff <- f
  }

  ## Check if g is specified as a function(x) only, then convert to a function(t,x)
  if(length(formals(g))==1){
    ggg <- function(t,x)g(x)
  }else
  {
    ggg <- g
  }

  ## Check if h is unspecified or specified as a function(x) only, then convert to a function(t,x)
  if(is.null(h)) h <- function(t,x)1
  if(length(formals(h))==1){
    hh <- function(t,x)h(x)
  }else
  {
    hh <- h
  }
  
  ## Check if r is unspecified or specified as a function(x) only, then convert to a function(t,x)
  if(is.null(r)) {
      r <- function(t,x)0
      NO.MORTALITY <- TRUE
  } else NO.MORTALITY <- FALSE
  
  if(length(formals(r))==1){
    rr <- function(t,x)r(x)
  }else
  {
    rr <- r
  }

  ## Find number of dimensions of the Brownian motion. Convert g, if necessary, to something that
  ## returns a nx-by-nB matrix
  g0 <- ggg(times[1],x0)

  if(is.null(dim(g0))) {
      nB <- 1
      gg <- function(t,x)matrix(ggg(t,x),nrow=nx,ncol=1)
  } else {
      nB <- ncol(g0)
      gg<-ggg
  }

  ## If the sample path of Brownian motion is not specified, then simulate it
  if(is.null(B)){
    B <- rvBM(times,nB)
  }

  if(!is.matrix(B)) B <- matrix(B,nrow=nt,ncol=nB)

  dB <- apply(B,2,diff)

  X <- array(NA,c(nt,nx))
  colnames(X) <- names(x0)
    X[1,] <- x0

  S <- numeric(nt)
  S[1] <- S0
  
  dt <- diff(times)

  for(i in 1:(nt-1))
  {
      ## Euler predictor
      fX <- as.numeric(ff(times[i],as.numeric(X[i,])))
      gX <- gg(times[i],X[i,])
      Y <- as.numeric(p(X[i,] + fX*dt[i] + as.numeric(gX %*% as.numeric(dB[i,]))))

      ## Corrector
      fY <- as.numeric(ff(times[i+1],Y))
      gY <- gg(times[i+1],Y)
      X[i+1,] <- as.numeric(p(X[i,] + 0.5*(fX+fY)*dt[i] + as.numeric(0.5*(gX+gY) %*% as.numeric(dB[i,]))))
      
      if(h(times[i+1],X[i+1,]) <= 0) break
      S[i+1] <- S[i] * exp(-rr(t,X[i,])*dt[i])
      if(S[i+1]<Stau) break
  }

  if(NO.MORTALITY)
  {
      return(list(times=times,X=X))
  }
  else {
      tau <- times[i+1]
      return(list(times=times[1:(i+1)],X=X[1:(i+1),],S=S[1:(i+1)],tau=tau))
  }      
}


#' Euler simulation of an Ito stochastic differential equation dX = f dt + g dB
#'
#' @name euler
#'
#' @param f function f(x) or f(t,x) giving the drift term in the SDE
#' @param g function g(x) or g(t,x) giving the noise term in the SDE
#' @param times numeric vector of time points, increasing
#' @param x0 numeric vector giving the initial condition
#' @param B sample path of (multivariate) Browniań motion. If omitted, a sample is drawn using rvBM
#' @param p Optional projection function that at each time poins projects the state, e.g. to enforce non-negativeness
#' @param h Optional stopping function; the simulation terminates if h(t,x)<=0
#' @examples
#' times <- seq(0,10,0.1)
#'
#' # Plot a sample path of the solution of the SDE dX=-X dt + dB
#' plot(times,euler(function(x)-x,function(x)1,times,0)$X,type="l")
#'
#' # Plot a sample path of the solution of the 2D time-varying SDE dX=(AX + FU) dt + G dB
#'
#' A <- array(c(0,-1,1,-0.1),c(2,2))
#' F <- G <- array(c(0,1),c(2,1))
#'
#' f <- function(t,x) A %*% x + F * sin(2*t)
#' g <- function(x) G
#'
#' BM <- rBM(times)
#'
#' plot(times,euler(f,g,times,c(0,0),BM)$X[,1],type="l")
#'
#' # Add solution corresponding to different initial condition, same noise
#'
#' lines(times,euler(f,g,times,c(1,0),BM)$X[,1],type="l",lty="dashed")
#' 
#' @export
euler <- function(f,g,times,x0,B=NULL,p=function(x)x,h=NULL,r=NULL,S0=1,Stau=runif(1))
{
  nx <- length(x0)
  nt <- length(times)

  ## Check if f is specified as a function(x) only, then convert to a function(t,x)
  if(length(formals(f))==1){
    ff <- function(t,x)f(x)
  }else
  {
    ff <- f
  }

  ## Check if g is specified as a function(x) only, then convert to a function(t,x)
  if(length(formals(g))==1){
    ggg <- function(t,x)g(x)
  }else
  {
    ggg <- g
  }

  ## Check if h is unspecified or specified as a function(x) only, then convert to a function(t,x)
  if(is.null(h)) h <- function(t,x)1
  if(length(formals(h))==1){
    hh <- function(t,x)h(x)
  }else
  {
    hh <- h
  }
  
  ## Check if r is unspecified or specified as a function(x) only, then convert to a function(t,x)
  if(is.null(r)) {
      r <- function(t,x)0
      NO.MORTALITY <- TRUE
  } else NO.MORTALITY <- FALSE
  
  if(length(formals(r))==1){
    rr <- function(t,x)r(x)
  }else
  {
    rr <- r
  }

  ## Find number of dimensions of the Brownian motion. Convert g, if necessary, to something that
  ## returns a nx-by-nB matrix
  g0 <- ggg(times[1],x0)

  if(is.null(dim(g0))) {
      nB <- 1
      gg <- function(t,x)matrix(ggg(t,x),nrow=nx,ncol=1)
  } else {
      nB <- ncol(g0)
      gg<-ggg
  }

  ## If the sample path of Brownian motion is not specified, then simulate it
  if(is.null(B)){
    B <- rvBM(times,nB)
  }

  if(!is.matrix(B)) B <- matrix(B,nrow=nt,ncol=nB)

  dB <- apply(B,2,diff)

  X <- array(NA,c(nt,nx))
  colnames(X) <- names(x0)
  X[1,] <- x0

  S <- numeric(nt)
  S[1] <- S0
  
  dt <- diff(times)

  for(i in 1:(nt-1))
  {
      fX <- as.numeric(ff(times[i],as.numeric(X[i,])))
      gX <- gg(times[i],as.numeric(X[i,]))
      
      X[i+1,] <- as.numeric(p( as.numeric(X[i,]) + fX*dt[i] + as.numeric(gX %*% as.numeric(dB[i,]))))

      if(h(times[i+1],X[i+1,]) <= 0) break
      S[i+1] <- S[i] * exp(-rr(t,X[i,])*dt[i])
      if(S[i+1]<Stau) break
  }

  colnames(X) <- names(x0)

  if(NO.MORTALITY)
  {
      return(list(times=times,X=X))
  }
  else {
      tau <- times[i+1]
      return(list(times=times[1:(i+1)],X=X[1:(i+1),],S=S[1:(i+1)],tau=tau))
  }      
}

#' Transition probabilities in a linear SDE dX = A*X*dt + u*dt + G*dB
#'
#'
#' @name dLinSDE
#'
#' @param A system matrix, quadratic numeric array or matrix
#' @param G noise matrix, same number of rows as A
#' @param t time, numeric scalar
#' @param x0 Initial state
#' @param u forcing. A numeric vector or matrix; see details
#' @param S0 Initial variance
#' @return A list containing St, the variance of the state at time t, and, depending on input arguments
#'
#'   eAt, the matrix that maps X0 to Xt, i.e. expm(A*t), provided x0==NULL
#'
#'   EX, the expectation of the state at time t, provided x0 is not NULL
#' @details
#' When returning EX, the expectation at time t, the input u can be specified as follows:
#'
#' * NULL, in which case it is assumed to be absent from the equation
#'
#' * a numeric vector with an element for each row in A. In this case u is assumed to be constant
#'
#' * a numeric array or matrix with same number of rows as A and two columns. In this case u is assumed
#'   to vary affinely in time, equaling the first column at time 0 and the second column at time t.
#'
#' The computations are not optimized for large systems, since they rely on the vector form of the matrix equations and do not use sparsity.
#' @examples
#' # A scalar equation with no input
#' (dLinSDE(-1,1,1))
#' # A scalar equation with constant input, starting in steady-state
#' (dLinSDE(A=-1,G=1,t=3,x0=1,u=1,S0=0.5))
#' 
#' @export
dLinSDE <- function(A,G,t,x0=NULL,u=NULL,S0=0*A)
{
    G <- as.matrix(G)
    A <- as.matrix(A)

    GG <- G %*% t(G)

    nx <- nrow(A)
    nB <- ncol(G)

    eAt <- as.matrix(Matrix::expm(A*t))

    St <- tryCatch(
        expr = {
            Sinf <- lyap(A,GG)
            St <- Sinf - eAt %*% (Sinf - S0) %*% t(eAt)            
        },
        error = function(e){
            I <- Matrix::Diagonal(n=nx,x=1)
            M <- A %x% I + I %x% A
            P <- rbind(array(0,c(nx^2,2*nx^2)),cbind(diag(rep(1,nx^2)),M))
            gs <- Matrix::expm(P*t) %*% c(as.numeric(GG),rep(0,nx^2))
            St <-eAt %*% S0 %*% t(eAt) + matrix(gs[nx^2+(1:(nx^2))],nrow=nx)
        })

    if(is.null(x0)){
      return(list(eAt = eAt,St = St))
    }

    if(is.null(u)){
      return(list(EX=eAt %*% x0,St=St))
    }

    ## Zero order hold
    if(length(u)==nx){
        I <- Matrix::Diagonal(n=nx,x=1)

        ## Extended state: (X,U)
        AA <- rbind(cbind(A,I),array(0,c(nx,2*nx)))
        Exu <- Matrix::expm(AA*t) %*% c(x0,u)
        EX <- Exu[1:nx]
        return(list(EX=EX,St=St))
    }

    ## If we reach this point, it must be first order hold
    u = as.matrix(u)

    ## First order hold
    if(all(dim(u)==c(nx,2))){
        O <- array(0,c(nx,nx))
        I <- Matrix::Diagonal(n=nx,x=1)

        ## Extended state: (dU/dt,U,X)
        AAA <- rbind(cbind(O,O,O),
                     cbind(I,O,O),
                     cbind(O,I,A))

        Exuu <- Matrix::expm(AAA*t) %*% c((u[,2]-u[,1])/t,u[,1],x0)
        EX <- tail(Exuu,nx)
        return(list(EX=EX,St=St))
    }
}

#' Algebraic Lyapunov equation A*X+X*t(A)+Q=0
#'
#' @name lyap
#'
#' @param A A quadratic matrix without eigenvalues on the imaginary axis
#' @param Q A symmetric matix of same dimension as A
#' @return X A symmetric matrix of same dimension as A
#'
#' @details
#' If A is asymptotically stable, Q is positive semidefinite and the pair (A,Q) is controllable, then X will be positive definite. Several similar results exist.
#' The implementation uses vectorization and kronecker products and does not employ sparsity, so is only suitable for small systems. 
#'
#' @examples
#' # A scalar example
#' (lyap(-1,1))
#' # A harmonic oscillator
#' (lyap(array(c(0,-1,1,-0.1),c(2,2)),diag(c(0,1))))
#'
#' @export
lyap <- function(A,Q) 
{
    A <- as.matrix(A)
    I <- diag(rep(1,nrow(A)))
    P <- kronecker(I,A)+kronecker(A,I)
    X <- -solve(P,as.numeric(Q))
    return(matrix(X,nrow=nrow(A)))
}

### Functionality for the Cox-Ingersoll-Ross process

## Internal function for converting CIR parameters to the parameters in the 
## non-central chi-squared distribution
cir.parameters <- function(x0,lambda,xi,gamma,t,Stratonovich=FALSE)
{
  ## If this is for the Stratonovich interpretation,
  ## convert first to the equivalent Ito equation
  if(Stratonovich) xi <- xi + gamma^2/4/lambda
  
  c <- 2*lambda/gamma^2/(1-exp(-lambda*t))
  nu <- 2*c*x0*exp(-lambda*t)
  n <- 4*lambda*xi/gamma^2
  
  return(list(c=c,nu=nu,n=n))
}

#' The Cox-Ingersoll-Ross process
#'
#' @description Density, distribution function, quantile function, and random generation 
#' for the transition probabilities in the Cox-Ingersoll-Ross process given by the 
#' stochastic differential equation dX = lambda*(xi-X)*dt + gamma*sqrt(X)*dB (interpreted 
#' in the sense of Ito (default) or Stratonovich (optional))
#'
#' @usage 
#' dCIR(x,x0,lambda,xi,gamma,t,Stratonovich=FALSE,log=FALSE)
#' 
#' @param x,q Target state, assumed >= 0
#' @param p Probability,  assumed >= 0 and <= 1.
#' @param x0 Initial state, assumed > 0
#' @param lambda Rate parameter, assumed > 0
#' @param xi Mean parameter, assumed > 0
#' @param gamma Noise intensity parameters, assumed > 0
#' @param t Terminal time, assumed > 0
#' @param Stratonovich Logical, TRUE for Stratonovich, FALSE (default) for Ito
#' @param log,log.p Logical, if TRUE, probabilities/densities are given as log(p). Default is FALSE
#' @param lower.tail Logical; if TRUE (default) probabilities are P(X<=x); otherise, P(X>x).
#'
#' @return dCIR gives the transition probability density, pCIR gives the distribution of the transitio probability, qCIR gives the quantiles, and rCIR samples a random terminal point.
#' 
#' The length of the result is determined by n for rCIR, and is the maximum of the lengths of the numerical arguments for the other functions. 
#'
#' @examples
#' x <- sort(rCIR(100,1,1,1,1,1))
#' par(mfrow=c(1,2))
#' plot(x,dCIR(x,1,1,1,1,1),ylab="p.d.f.")
#' F <- pCIR(x,1,1,1,1,1)
#' plot(x,F)
#' lines(qCIR(F,1,1,1,1,1),F)
#' @export
dCIR <- function(x,x0,lambda,xi,gamma,t,Stratonovich=FALSE,log=FALSE)
{
  p <- cir.parameters(x0,lambda,xi,gamma,t,Stratonovich)
  
  ## The following is the explicit expression
  ## d <- c*exp(-0.5*(2*c*x+nu))*(2*c*x/nu)^(n/4-0.5)*besselI(sqrt(2*c*nu*x),n/2-1)
  ## More accurately, and easier, to use the  built-in p.d.f. of the
  ## non-central chi-squared distribution
  ld <- log(2*p$c) + dchisq(2*p$c*x,p$n,p$nu,log=TRUE)
  
  if(log==TRUE) return(ld) else return(exp(ld))
}

#' @rdname dCIR
#' @usage 
#' pCIR(x,x0,lambda,xi,gamma,t,Stratonovich=FALSE,log.p=FALSE,lower.tail=TRUE)
#' @export
pCIR <- function(x,x0,lambda,xi,gamma,t,Stratonovich=FALSE,lower.tail=TRUE,log.p=FALSE)
{
  p <- cir.parameters(x0,lambda,xi,gamma,t,Stratonovich)
  
  lp <- pchisq(2*p$c*x,p$n,p$nu,lower.tail=lower.tail,log.p=log.p)
  
  return(lp)
}

#' @rdname dCIR
#' @usage
#' qCIR(p,x0,lambda,xi,gamma,t,Stratonovich=FALSE,log.p=FALSE,lower.tail=TRUE)
#' @export
qCIR <- function(ps,x0,lambda,xi,gamma,t,Stratonovich=FALSE,lower.tail=TRUE,log.p=FALSE)
{
  p <- cir.parameters(x0,lambda,xi,gamma,t,Stratonovich)
  
  x <- qchisq(ps,p$n,p$nu,lower.tail=lower.tail,log.p=log.p)/2/p$c
  
  return(x)
}

#' @rdname dCIR
#' @usage
#' rCIR(n,x0,lambda,xi,gamma,t,Stratonovich=FALSE)
#' @export
rCIR <- function(n,x0,lambda,xi,gamma,t,Stratonovich=FALSE)
{
  p <- cir.parameters(x0,lambda,xi,gamma,t,Stratonovich)
  
  x <- rchisq(n,p$n,p$nu)/2/p$c
  
  return(x)
}

## Internal function for converting OU parameters to the parameters in the 
## Gaussian distribution
ou.parameters <- function(x0,lambda,xi,sigma,t)
{
  mean <- xi + (x0-xi)*exp(-lambda*t)
  sd <- if(lambda==0) sd <- sigma*sqrt(t) else sd <- sigma*sqrt((1-exp(-2*lambda*t))/2/lambda)
  
  return(list(mean=mean,sd=sd))
  }

#' The Ornstein-Uhlenbeck process
#'
#' @description Density, distribution function, quantile function, and random generation 
#' for the transition probabilities in the (shifted) Ornstein-Uhlenbeck process given by the 
#' stochastic differential equation dX = lambda*(xi-X)*dt + sigma*dB 
#'
#' @usage 
#' dOU(x,x0,lambda,xi,sigma,t,log=FALSE)
#' 
#' @param x,q Target state, assumed >= 0
#' @param p Probability,  assumed >= 0 and <= 1.
#' @param x0 Initial state, assumed > 0
#' @param lambda Rate parameter, assumed > 0
#' @param xi Mean parameter, assumed > 0
#' @param sigma Noise intensity parameters, assumed > 0
#' @param t Terminal time, assumed > 0
#' #' @param log,log.p Logical, if TRUE, probabilities/densities are given as log(p). Default is FALSE
#' @param lower.tail Logical; if TRUE (default) probabilities are P(X<=x); otherise, P(X>x).
#'
#' @return dOU gives the transition probability density, pOU gives the distribution of the transitio probability, qOU gives the quantiles, and rOU samples a random terminal point.
#' 
#' The length of the result is determined by n for rOU, and is the maximum of the lengths of the numerical arguments for the other functions. 
#'
#' @examples
#' x <- sort(rOU(100,1,1,1,1,1))
#' par(mfrow=c(1,2))
#' plot(x,dOU(x,1,1,1,1,1),ylab="p.d.f.")
#' F <- pOU(x,1,1,1,1,1)
#' plot(x,F)
#' lines(qOU(F,1,1,1,1,1),F)
#' @export
dOU <- function(x,x0,lambda,xi,sigma,t,log=FALSE)
{
  p <- ou.parameters(x0,lambda,xi,sigma,t)
  
  ld <- dnorm(x,p$mean,p$sd,log=TRUE)
  
  if(log==TRUE) return(ld) else return(exp(ld))
}

#' @rdname dOU
#' @usage 
#' pOU(x,x0,lambda,xi,sigma,t,log.p=FALSE,lower.tail=TRUE)
#' @export
pOU <- function(x,x0,lambda,xi,sigma,t,lower.tail=TRUE,log.p=FALSE)
{
  p <- ou.parameters(x0,lambda,xi,sigma,t)
  
  lp <- pnorm(x,p$mean,p$sd,lower.tail=lower.tail,log.p=log.p)
  
  return(lp)
}

#' @rdname dOU
#' @usage
#' qOU(p,x0,lambda,xi,sigma,t,log.p=FALSE,lower.tail=TRUE)
#' @export
qOU <- function(ps,x0,lambda,xi,sigma,t,lower.tail=TRUE,log.p=FALSE)
{
  p <- ou.parameters(x0,lambda,xi,sigma,t)
  
  x <- qnorm(ps,p$mean,p$sd,lower.tail=lower.tail,log.p=log.p)
  
  return(x)
}

#' @rdname dOU
#' @usage
#' rOU(n,x0,lambda,xi,sigma,t)
#' @export
rOU <- function(n,x0,lambda,xi,sigma,t)
{
  p <- ou.parameters(x0,lambda,xi,sigma,t)
  
  x <- rnorm(n,p$mean,p$sd)
  
  return(x)
}

##
