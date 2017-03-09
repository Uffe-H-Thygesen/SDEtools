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

