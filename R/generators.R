#' Discretize scalar advection-diffusion equation using finite volumes
#'
#' @description
#' fvade discretizes the advection-diffusion equation
#'
#'      dC/dt = -(u*C-D*C')'
#'
#'   on an interval [a,b] using the finite-volume method.
#' @param  u function mapping state (numeric scalar) to advective term (numeric scalar)
#' @param  D function mapping state (numeric scalar) to diffusivity (numeric scalar)
#' @param  xgrid The numerical grid. Numeric vector of increasing values, giving cell boundaries
#' @param  bc String indicating boundary conditions. See details.
#' @return a quadratic matrix, the generator of the approximating continuous-time Markov chain, with length(xgrid)-1 columns
#'
#' @details Handling of boundary conditions: Input argument bc is a single character coding boundary conditions as follows:
#'          'r': Reflecting boundaries
#'          'p': Periodic boundaries
#'          'a': Absorbing boundaries. In this case G will be a sub-generator
#'          'e': Return generator, Extended to include absorbing boundaries
#'
#'  Return value: The function fvade returns a generator (or sub-generator) G of a continuous-time Markov Chain. This chain jumps
#'  between cells defined by xgrid. When using the generator to solve the Kolmogorov equations, note that G operates on
#'  probabilities of each cell, not on the probability density in each cell. The distinction is particularly important when the
#'  grid is non-uniform.
#'
#' @examples
#' # Generator of standard Brownian motion with unit drift on the unit interval, reflecting boundaries
#' xi <- seq(0,1,0.05)    # Cell boundaries
#' dx <- diff(xi)         # Cell widths
#' xc <- xi[-1] - 0.5*dx  # Cell centers
#' G <- fvade(function(x)1,function(x)0.5,seq(0,1,0.05),'r')
#'
#' # Find the density of the stationary distribution
#' phi <- Null(G)         # Find unnormalized stationary probabilities
#' phi <- phi/sum(phi)/dx # Convert to densities
#' plot(xc,phi,type="l",xlab="x",ylab="Stationary density")
#'
#' @export
fvade <- function(u,D,xgrid,bc)
{
  # Number of grid cells
  nc = length(xgrid)-1;
  dx = diff(xgrid);

  # Indeces of left and right cell, with wrap around
  Ileft = c(nc,1:(nc-1))
  Iright= c(2:nc,1)

  # Diffusivities at left interfaces
  Dil = sapply(xgrid[1:nc],D);
  Dir = sapply(xgrid[2:(nc+1)],D);

  # Diffusion elements of the generator
  Dl = 2*Dil/dx/(dx + dx[Ileft]);
  Dr = 2*Dir/dx/(dx + dx[Iright]);

  # Advection at left and right interfaces
  Uil = sapply(xgrid[1:nc],u);
  Uir = sapply(xgrid[2:(nc+1)],u);

  # Advection elements of the generator
  Ul = pmax(-Uil,0)/dx;
  Ur = pmax( Uir,0)/dx;

  # Mimic Matlab for off-diagonals
  mydiag <- function(v,offset)
  {
    Dv <- as.array(diag(v))
    n <- length(v)
    if(offset==0) return(Dv)
    if(offset< 0) return(cbind(rbind(array(0,c(-offset,n)),Dv),array(0,c(n-offset,-offset))))
    if(offset> 0) return(cbind(array(0,c(n+offset,offset)),rbind(Dv,array(0,c(offset,n)))))
  }

  # Extended generator
  Ge = mydiag(c(Dl+Ul,0),-1) + mydiag(c(0,Dr+Ur),1);

  # Handle boundary conditions

  if(bc=='p')
  {
    G = Ge[2:(nc+1),2:(nc+1)]
    G[1,nc] = Ge[2,1]
    G[nc,1] = Ge[nc-1,nc]

    G = G - diag(apply(G,1,sum))
  }


  if(bc=='r')
  {
    G = Ge[2:(nc+1),2:(nc+1)];

    G = G - diag(apply(G,1,sum))
  }

  if(bc=='a')
  {
    Ge = Ge - diag(apply(Ge,1,sum))

    G = Ge[2:(nc+1),2:(nc+1)];
  }

  if(bc=='e')
  {
    Ge = Ge - diag(apply(Ge,1,sum))
    G = Ge;
  }

  return(G)
}
