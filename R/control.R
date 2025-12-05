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

#' Solve the optimal control problem using policy iteration, i.e. find the optimal strategy
#' and value function for the case where the uncontrolled system is given by a subgenerator G0
#' (either due to discounting or due to absorbing boundaries)
#'
#' @param  G0 The sub-generator of the uncontrolled system 
#' @param  G1 A list of (sub-)generators for each control
#' @param  k The running cost 
#' @param  uopt A list of functions returning optional controls in each state. See details. 
#' @param  iter.max = 1000 Maximum number of iterations
#' @param  tol = 1e-12 Tolerance for convergence 
#' @param  do.minimize=TRUE Is the optimization problem one of minimization (and not maximization)?
#' @param  do.return.QSD=FALSE Compute and return the quasi-stationary distribution 
#' @return A list containing
#' 
#' V: The value function, as a vector with an element for each state
#' 
#' u: The optimal controls, as a matrix with a row for each state and a column for each control
#'
#' If do.return.QSD==TRUE, also:
#' 
#' qsd.value: The decay rate of the quasi-stationary distribution
#' 
#' qsd.vector: The quasi-stationary distribution
#'
#' @details Specifying optimal controls: Each control variable has an element in the list G1
#' of (sub)generators, and an element in the list uopt of optimal controls.
#'
#' The function uopt[[i]] is applied to G1[[i]]*V, where V is the value function. The result
#' is the optimal value of the control, in each possible state.
#'
#' The final generator for the controlled system is obtain as G0 + sum(Diagonal(u[[i]])*G1[[i]])
#' where the sum is over the controls i.
#' 
#' @examples
#' ## Controlling a system to the boundary with minimum effort
#' H <- 2
#' xi <- seq(-H,H,length=201)
#' xc <- as.numeric(cell.centers(xi,c(0,1))$x)
#' dx <- diff(xi)
#' 
#' G0 <- fvade(function(x)0,function(x)0.5,xi,'a')
#' Gp <- fvade(function(x)1,function(x)0,xi,'a')   # Going right 
#' Gn <- fvade(function(x)-1,function(x)0,xi,'a')  # Going left
#'
#' # The payoff - time plus control effort
#' k <- function(u) 0.5 + 0.5*u[,1]^2 + 0.5*u[,2]^2
#'
#' # For each control, the optimal u maximizes Li V*u + k over [0,Inf)
#' # Here, going right corresponds to Lp V = V' while going left corresponds to
#' # Ln V = -V'
#' uopt <- function(dV)pmax(0,-dV)
#'
#' sol <- PolicyIterationRegular(G0,list(Gp,Gn),k,list(uopt,uopt),do.return.QSD=TRUE)
#' 
#' par(mfrow=c(1,3))
#' plot(xc,sol$V,xlab="x",ylab="Value function")
#' lines(xc,log(cosh(H)/cosh(xc)))
#' plot(xc,sol$u[,1]-sol$u[,2],xlab="x",ylab="Optimal control")
#' lines(xc,tanh(xc))
#' plot(xc,sol$qsd.vector/dx,type="l",xlab="x",ylab="QSD",main=sol$qsd.value)
#' @export
PolicyIterationRegular <- function(G0,G1,k,uopt,iter.max = 1000,tol = 1e-12,
                                   do.minimize=TRUE,do.return.QSD=FALSE)
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
        rhs <- -k(u)
        V <- as.numeric(Matrix::solve(Gcl,rhs))

        ## Determine optimal strategy in response to the current value
        for(i in 1:nu) u[,i] <- uopt[[i]](as.numeric(G1[[i]] %*% V))
        
        if( xor(max(V-Vold)<tol,do.minimize)) break
        if(iter > iter.max) break
        
        Vold <- V
        iter <- iter + 1
    }

    if(do.return.QSD)
    {
        qsd <- QuasiStationaryDistribution(Gcl)
        return(list(V=V,u=u,
                    qsd.value=qsd$value,
                    qsd.vector=qsd$vector))
    }
    
    return(list(V=V,u=u))
}

#' Solve the optimal control problem using policy iteration, i.e. find the optimal strategy
#' and value function for the case where the uncontrolled system is given by a generator G0
#'
#' @param  G0 The generator of the uncontrolled system 
#' @param  G1 A list of generators for each control
#' @param  k The running cost 
#' @param  uopt A list of functions returning optional controls. See details. 
#' @param  iter.max = 1000 Maximum number of iterations
#' @param  tol = 1e-12 Tolerance for convergence 
#' @param  do.minimize=TRUE
#' @return A list containing
#' V: The value function, as a vector with an element for each state
#' u: The optimal controls, as a matrix with a row for each state and a column for each control
#' pi: The stationary distribution for the controlled system
#' gamma: The expected running cost in stationarity
#'
#' @details Specifying optimal controls: Each control variable has an element in the list G1
#' of generators, and an element in the list uopt of optimal controls.
#'
#' The function uopt[[i]] is applied to G1[[i]]*V, where V is the value function. The result
#' is the optimal value of the control, in each possible state.
#'
#' The final generator for the controlled system is obtain as G0 + sum(Diagonal(u[[i]])*G1[[i]])
#' where the sum is over the controls i.
#'
#' @examples
#' require(SDEtools)
#'
#' ## Fisheries management dX = X*(1-X)*dt - U*dt + sigma*X*dB
#' sigma <- 0.5
#' ## Rewrite in advection-diffusion form
#' D <- function(x) 0.5*sigma^2*x^2
#' Dp <- function(x) sigma^2*x
#' u <- function(x) x*(1-x) - Dp(x)
#' 
#' # Spatial grid - concentrated around 0
#' xi <- seq(0,sqrt(3),length=201)^2
#' xc <- 0.5*(head(xi,-1) + tail(xi,-1))
#' 
#' # Generator for the unfished system, and for fishing
#' G0 <- fvade(u,D,xi,'r')
#' G1 <- fvade(function(x)-1,function(x)0,xi,'r')
#'
#' # Pay-off function
#' k <- function(u) sqrt(u)
#' 
#' # Optimal control - really 1/4/V'(x)^2, but avoid V'<=0, and do not
#' # fish in the lowest element in state space, to avoid degeneracies
#' uopt <- function(dV) { u <- 0.25/pmax(1e-6,-dV)^2; u[1] <- 0 ; return(u)}
#'
#' sol <- PolicyIterationSingular(G0,list(G1),k,list(uopt),do.minimize=FALSE)
#'
#' # Plots comparing with analytical solution: V(x) = 0.5*log(x), u(x)=x^2
#' par(mfrow=c(1,2))
#' i1 <- sum(xi<1) # Identify state corresponding to x=1, for shifting V
#' plot(xc,0.5*log(xc),type="l",ylim=c(-2,0.6),xlab="x",ylab="Value function")
#' points(xc,sol$V - sol$V[i1])
#' plot(xc,xc^2,type="l",xlab="x",ylab="Control")
#' points(xc,sol$u)
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
    V <- numeric(nrow(G0))

    e <- rep(1,nrow(G0))

    iter <- 1
    
    while(TRUE)
    {
        ## Determine optimal strategy in response to the current value
        for(i in 1:nu) u[,i] <- uopt[[i]](as.numeric(G1[[i]] %*% V))

        ## Construct closed loop generator
        Gcl <- G0
        for(i in 1:length(G1)) Gcl <- Gcl + Diagonal(x=u[,i]) %*% G1[[i]]

        ## Extend with ones right and below; have 0 in the bottom right
        Ge <- rbind(cbind(Gcl,-e),c(e,0))

        ## Solve for the value function and the average performance
        rhs <- c(-k(u),0)
        Vg <- Matrix::solve(Ge,rhs)
        V <- as.numeric(head(Vg,-1))
        gamma <- as.numeric(tail(Vg,1))

        if(xor(gamma-gammaold<tol,do.minimize)) break
        if(iter > iter.max) break

        gammaold <- gamma
        iter <- iter + 1
    }

    pi <- StationaryDistribution(Gcl)
    
    return(list(V=V,u=u,pi=pi,gamma=gamma))
}
