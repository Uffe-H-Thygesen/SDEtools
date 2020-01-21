Riccati <- function(time,Ss,par)
{
  S <- matrix(Ss[1:(length(Ss)-1)],nrow=nrow(par$A))
  SA <- S %*% par$A
  dSdt <- SA + t(SA) - S %*% par$FRiF %*% S + par$Q
  dsdt <- 0.5*sum(diag(par$Gt %*% S %*% par$G))

  return(list(c(as.numeric(dSdt),dsdt)))
}

#' Solve the time-varying LQR (Linear Quadratic Regulator) problem
#'
#' @param  times numeric vector of time points
#' @param  A system matrix, n-by-n numeric array
#' @param  F control matrix, n-by-m numeric array
#' @param  G noise matrix, n-by-l numeric array
#' @param  Q state penalty running cost, n-by-n numeric array specifying the quadratic form in x in the running cost
#' @param  R control penalty in running cost, m-by-m numeric array specifying the quadratic form in u in the running cost
#' @param  P terminal state penalty, n-by-n numeric array specifying the quadratic form in x in the terminal cost
#' @return A list containing
#' times, A, F, G, Q, R, P: The input arguments
#' S the value function. A length(times)*n*n array containing, for each time point, the quadratic form in x in the value function
#' s the value function. A length(times) vector containing, for each time point, the off-set in the value function, i.e. the value for x=0
#' L the optimal control. A length(times)*m*n array containing, for each time point, the gain matrix Lt that maps states to optimal controls
#'
#' @examples
#' # Feedback control of a harmonic oscillator
#' A <- array(c(0,1,-1,0),c(2,2))
#' F <- G <- array(c(0,1),c(2,1))
#' Q <- diag(rep(1,2))
#' R <- 1
#' times <- seq(0,10,0.1)
#' sol <- lqr(times,A,F,G,Q,R)
#' matplot(-times,array(sol$S,c(length(times),length(A))),type="l",xlab="Time",ylab="Elements in S")
#' legend("topright",c("xx","xv","vx","vv"),lty=c("solid","dashed","dashed","dotted"),col=1:4)
#' matplot(-times,array(sol$L,c(length(times),length(F))),type="l",
#'          xlab="Time",ylab="Feedback control",lty=c("solid","dashed"))
#' legend("topright",c("Position","Velocity"),lty=c("solid","dashed"),col=1:2)
#' @export
lqr <- function(times,A,F,G,Q,R,P=NULL)
{
  RiF <- solve(R,t(F))

  if(is.null(P)) P <- array(0,rep(nrow(A),2))

  par <- list(A = A,G = G,Gt = t(G),FRiF = F %*% RiF,Q=Q)

  require(deSolve)

  sol <- ode(y = c(as.numeric(P),0),times=times,func=Riccati,parms=par)

  SvecToGain <- function(Ss) as.numeric(-RiF %*% matrix(Ss[2:(length(Ss)-1)],nrow=nrow(par$A)))

  s <- sol[,ncol(sol)]
  S <- array(sol[,2:(ncol(sol)-1)],c(length(times),nrow(A),nrow(A)))
  L <- apply(sol,1,SvecToGain)

  L <- array(t(L),c(length(times),ncol(F),nrow(A)))

  return(list(times=times,A=A,F=F,G=G,Q=Q,R=R,P=P,S=S,s=s,L=L))
}

## Find the optimal strategy and value function for the case where
## the uncontrolled system is given by a subgenerator G0
## (either due to discounting or due to absorbing boundaries)
#' @export
PolicyIterationRegular <- function(G0,G1,k,uopt,iter.max = 1000,tol = 1e-12,do.minimize=TRUE)
{
    if(!is.list(G1)) G1 <- list(G1)
    if(!is.list(uopt)) uopt <- list(uopt)

    Vold <- (2*do.minimize-1)*Inf

    ## Initial guess on u obtained with a zero value function 
    nu <- length(uopt)
    u <- array(NA,c(nrow(G0),nu))
    for(i in 1:nu) u[,i] <- uopt[[i]](numeric(nrow(G0)))

    iter <- 1
    
    while(TRUE)
    {
        ## Evaluate performance of current strategy
        Gcl <- G0
        for(i in 1:length(G1)) Gcl <- Gcl + Diagonal(x=u[,i]) %*% G1[[i]]
        V <- as.numeric(solve(Gcl,-k(u)))

        ## Determine optimal strategy in response to the current value
        for(i in 1:nu) u[,i] <- uopt[[i]](as.numeric(G1[[i]] %*% V))
        
        if( xor(print(max(V-Vold))<tol,do.minimize)) break
        if(iter > iter.max) break
        
        Vold <- V
        iter <- iter + 1
    }

    return(list(V=V,u=u))
}

## Find the optimal strategy and value function for the case where
## the uncontrolled system is given by a generator G0
#' @export
PolicyIterationSingular <- function(G0,G1,k,uopt,iter.max = 1000,tol = 1e-12,do.minimize=TRUE)
{
    if(!is.list(G1)) G1 <- list(G1)
    if(!is.list(uopt)) uopt <- list(uopt)

    ## Initial guess on the optimal performance
    gammaold <- (2*do.minimize-1)*Inf

    ## Initial guess on u obtained with a zero value function 
    nu <- length(uopt)
    u <- array(NA,c(nrow(G0),nu))
    for(i in 1:nu) u[,i] <- uopt[[i]](numeric(nrow(G0)))

    e <- rep(1,nrow(G0))

    iter <- 1
    
    while(TRUE)
    {
        ## Construct closed loop generator
        Gcl <- G0
        for(i in 1:length(G1)) Gcl <- Gcl + Diagonal(x=u[,i]) %*% G1[[i]]

        ## Extend with ones right and below; have 0 in the bottom right
        Ge <- rbind(cbind(Gcl,-e),c(e,0))

        ## Solve for the value function and the average performance
        Vg <- solve(Ge,c(-k(u),0))
        V <- as.numeric(head(Vg,-1))
        gamma <- as.numeric(tail(Vg,1))

        ## Determine optimal strategy in response to the current value
        for(i in 1:nu) u[,i] <- uopt[[i]](as.numeric(G1[[i]] %*% V))

        if(xor(print(max(gamma-gammaold))<tol,do.minimize)) break
        if(iter > iter.max) break

        gammaold <- gamma
        iter <- iter + 1
    }

    return(list(V=V,u=u,gamma=gamma))
}
