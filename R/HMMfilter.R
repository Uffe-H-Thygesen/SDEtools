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
##' @param x0dist Prior distribution of the initial state. See details. 
##' @param tvec Vector of (increasing) time points where observations are taken
##' @param yvec Vector of observations at each time point
##' @param lfun Likelihood function so that lfun(x,y) gives the likelihood of y given x
##' @param do.smooth Do we want smoothing, or only predictive filtering / estimation?
##' @param N.sample Number of "typical tracks" sampled (defaults to 0)
##' @param do.Viterbi Do we want the most probable state sequence, found with the Viterbi algorithm?
##' @param pfun C.d.f. of observations given states, i.e. pfun(x,y) gives P(Y<=y | X = x). If supplied, pseudo-prediction residuals will be computed
##' 
##' @details 
##' The distribution of the initial condition x0 can be specified in a number of ways: 
##'   . If x0dist is a function, it is interpreted as the c.d.f. of the initial state
##'   . If x0dist is a single number, it is interpreted as a deterministic initial state
##'   . If x0dist is a numeric vector of the same length as xi, then it is interpreted as a c.d.f.
##'   . If x0dist is a numeric vector with one less element than xi, then it is interpreted as cell probabilities.
##'   
##' @return A list containing:
##' phi A tabulation of the predicitive probability densities
##' psi A tabulation of the estimated probability densities
##' pi  (If do.smooth==TRUE) A tabulation of the smoothed probability densities
##' c   A vector containing the normalization constants at each time step
##' loglik The total log-likelihood of the model
##' Xmpt A vector containing the most probabile path (if do.Viterbi==TRUE)
##' U A vector containint the pseudo-prediction residuals (if pfun is supplied)
##' 
##' @author Uffe Høgsbro Thygesen
##'
##' @export
HMMfilterSDE <-
    function(u,D,xi,bc,x0dist,tvec,yvec,lfun,do.smooth=FALSE,N.sample=0,do.Viterbi=FALSE,pfun=NULL)
{
    G <- fvade(u,D,xi,bc);

    xc <- cell.centers(xi)

    nx <- length(xc)
    nt <- length(tvec)

    dt <- diff(tvec)

    if(var(dt)/(mean(dt)^2) < 1e-8)
    {
        dt <- mean(dt)
    } else {
        stop("Non-uniform sampling has not been implemented yet.")
    }

    ltab <- outer(xc,yvec,lfun)

    phi <- array(0,c(nt,nx))
    psi <- phi
    c <- numeric(nt)
    
    ## Set initial distribution
    phi0 <- NULL
    if(is.function(x0dist))
    {
        phi0 <- diff(x0dist(xi));
        phi0 <- phi0 / sum(phi0)
    }
    if(is.numeric(x0dist))
    {
        if(length(x0dist)==1) phi0 <- diff(xi>x0dist)       # x0dist contains a deterministic IC
        if(length(x0dist)==nx) phi0 <- x0dist               # x0dist contains probabilities of cells
        if(length(x0dist)==length(xi)) phi0 <- diff(x0dist) # x0dist contains c.d.f. at cell interfaces
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

    for(i in 1:(nt-1))
    {
        psi[i,] = phi[i,] * ltab[,i];      # Data update
        c[i] <- sum(psi[i,])
        psi[i,] = psi[i,] / c[i]           # Normalize
        phi[i+1,] = psi[i,] %*% P;         # Time update
    }
    i <- i + 1
    psi[i,] = phi[i,] * ltab[,i];      # Data update
    c[i] <- sum(psi[i,])
    psi[i,] = psi[i,] / c[i]           # Normalize
    loglik = sum(log(c))
    
    ans <- list(phi=phi,psi=psi,c=c,loglik=loglik)

    if(do.smooth) {
        mu <- array(NA,c(length(xc),length(tvec)))
        lambda <- mu 

        pi <- t(mu)

        lambda[,length(tvec)] <- ltab[,length(tvec)]
        cvec <- rep(NA,length(tvec))

        for(i in length(tvec):2)
        {
            mu[,i-1] <- as.numeric(P %*% lambda[,i])
            lambda[,i-1] <- mu[,i-1] * ltab[,i-1];

            cvec[i-1] <- max(lambda[,i-1])
  
            lambda[,i-1] <- lambda[,i-1] / cvec[i-1]
        }

        for(i in 1:length(tm))
        {
            pi[i,] <- phi[i,] * lambda[,i]
            pi[i,] <- pi[i,] / sum(pi[i,])
        }

        ans <- c(ans,list(pi=pi,lambda=lambda,mu=mu))
    }

    if(do.Viterbi) {
        ## Backward dynamic programming iteration
        ## V is the value function (log-likelihood);
        ## U is the strategy (next state)
        U <- V <- array(NA,c(nx,nt))

        lP <- log(P)
        
        V[,nt] <- log(ltab[,nt])
        for(t in nt:2)
        {
            W <- outer(log(ltab[,t-1]),rep(1,nx)) + lP + outer(rep(1,nx),V[,t])
            V[,t-1] <- apply(W,1,max)
            U[,t-1] <- apply(W,1,which.max)
        }

        ## Forward loop, finding the path
        X <- numeric(nt)
        
        X[1] <- which.max(V[,1])
        for(t in 2:nt) X[t] <- U[X[t-1],t-1]

        ans <- c(ans,list(Xmpt=xc[X]))
    }


    if(N.sample>0)
    {
        Xtyp <- array(NA,c(nt,N.sample))

        for(j in 1:N.sample)
        {
            Xtyp[nt,j] <- sample.int(nx,size=1,prob=psi[nt,])
            for(t in nt:2)
            {
                prob <- as.numeric(psi[t-1,]) * as.numeric(P[,Xtyp[t,j]])
                prob <- prob / sum(prob)
                Xtyp[t-1,j] <- sample.int(nx,size=1,prob=prob)
            }
        }
        
        Xtyp <- apply(Xtyp,2,function(x)xc[x])

        ans <- c(ans,list(Xtyp = Xtyp))
    }
    
    if(!is.null(pfun))
    {
        U <- numeric(length(tvec))
        for(i in 1:length(tvec))
        {
            U[i] <- sum(phi[i,]*pfun(xc,yvec[i]))
        }

        ans <- c(ans,list(U=U))
    }

    return(ans)
}
