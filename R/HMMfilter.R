##' Estimate states in a scalar stochastic differential equation based on discretization to a HMM
##'
##' HMMfilterSDE computes the state estimation for a discretely observed stochastic differential equation
##'
##' @name HMMfilterSDE
##' 
##' @param u function mapping state (numeric scalar) to advective term (numeric scalar)
##' @param D function mapping state (numeric scalar) to diffusivity (numeric scalar)
##' @param xi The numerical grid. Numeric vector of increasing values, giving cell boundaries
##' @param bc String indicating boundary conditions. See details.
##' @param x0dist Prior distribution of the initial state, specified as a c.d.f. or as a vector of same length as xi, start
##' @param tvec Vector of (increasing) time points where observations are taken
##' @param yvec Vector of observations at each time point
##' @param lfun Likelihood function so that lfun(x,y) gives the likelihood of y given x
##' @param do.smooth Do want smoothing, or only predictive filtering / estimation?
##' @param do.Viterbi Do we want the most probable state sequence, found with the Viterbi algorithm?
##' @param pfun C.d.f. of observations given states, i.e. pfun(x,y) gives P(Y<=y | X = x). If supplied, pseudo-prediction residuals will be computed
##' 
##' @return a quadratic matrix, the generator of the approximating continuous-time Markov chain, with length(xgrid)-1 columns
##' 
##' @author Uffe Høgsbro Thygesen
##'
##' @export
HMMfilterSDE <-
    function(u,D,xi,bc,x0dist,tvec,yvec,lfun,do.smooth=FALSE,do.Viterbi=FALSE,pfun=NULL)
{
    G <- fvade(u,D,xi,bc);

    xc <- head(xi,-1) + 0.5*diff(xi)

    nx <- length(xc)
    nt <- length(tvec)

    dt <- diff(tvec)

    if(var(dt)/(mean(dt)^2) < 0.01)
    {
        dt <- mean(dt)
    } else {
        stop("Non-uniform sampling has not been implemented yet.")
    }

    ltab <- outer(xc,yvec,lfun)

    phi <- array(0,c(nt,nx))
    psi <- phi

    ## Set initial distribution
    phi0 <- NULL
    if(is.function(x0dist))
    {
        phi0 <- diff(x0dist(xi));
        phi0 <- phi0 / sum(phi0)
    }
    if(is.numeric(x0dist))
    {
        if(length(x0dist)==length(xc)) phi0 <- x0dist        # x0dist contains probabilities of cells
        if(length(x0dist)==length(xi)) phi0 <- diff(x0dist)  # x0dist contains c.d.f. at cell interfaces
    }
        
    if(is.null(x0dist))
    {
        phi0 <- StationaryDistribution(G)
    }

    if(is.null(phi0)) stop("Initial distribution not specified correctly")
    if(any(phi0<0)) stop("Some initial probability is negative")
    if(abs(sum(phi0)-1)>1e-7) stop("Initial probabilities do not sum to one")

    phi[1,] <- phi0
    
    ## Transition probabilities 
    P <- as.matrix(Matrix::expm(G*dt))

    for(i in 1:(length(tvec)-1))
    {
        psi[i,] = phi[i,] * ltab[,i];      # Data update
        psi[i,] = psi[i,] / sum(psi[i,])   # Normalize
        phi[i+1,] = psi[i,] %*% P;         # Time update
    }
    i <- i + 1
    psi[i,] = phi[i,] * ltab[,i];      # Data update
    psi[i,] = psi[i,] / sum(psi[i,])   # Normalize

    
    if(do.smooth) print("Smoothing has not been implemented yet.")

    if(do.Viterbi) print("The Viterbi algorithm has not been implemented yet.")

    if(!is.null(pfun)) print("Psuedo-residuals have not been implemented yet.")

    return(list(phi=phi,psi=psi))
}
