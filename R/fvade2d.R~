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
#' @param  sparse logical indicating if the result should be returned as a sparse matrix
#' @return a quadratic matrix, the generator of the approximating continuous-time Markov chain, with length(xgrid)-1 columns
#'
#' @details Handling of boundary conditions: Input argument bc is a single character, or a vector of two characters, coding the condition at each boundary as follows:
#'          'r': Reflecting boundary
#'          'p': Periodic boundaries: Exit at this boundary to re-enter at the other
#'          'a': Absorbing boundaries. In this case G will be a sub-generator
#'          'c': Continue beyond boundary (experimental; read the source)
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
fvade <- function(u,D,xgrid,bc,sparse=TRUE)
{
    require(Matrix)

    if(length(bc)==1) bc <- rep(bc,2)
    
    ## Number of grid cells
    nc = length(xgrid)-1;
    dx = diff(xgrid);

    ## Indeces of left and right cell, with wrap around
    Ileft = c(nc,1:(nc-1))
    Iright= c(2:nc,1)

    ## Diffusivities at interfaces
    Di = sapply(xgrid,D);

    ## Diffusion elements of the generator
    Dl = 2*head(Di,-1)/dx/(dx + dx[Ileft]);
    Dr = 2*tail(Di,-1)/dx/(dx + dx[Iright]);

    ## Advection at left and right interfaces
    Ui <- sapply(xgrid,u);
    Uil = head(Ui,-1)
    Uir = tail(Ui,-1)

    ## Advection elements of the generator
    Ul = pmax(-Uil,0)/dx;
    Ur = pmax( Uir,0)/dx;

    ## Mimic Matlab for off-diagonals
    mydiag <- function(v,offset=0)
    {
        size <- length(v)+abs(offset)
        j <- (1:length(v)) + max(0,offset)
        i <- (1:length(v)) + max(0,-offset)
        return(sparseMatrix(i=i,j=j,x=v,dims=rep(size,2)))
    }
    
    Ge = mydiag(c(Dl+Ul,0),-1) + mydiag(c(0,Dr+Ur),1) -mydiag(c(0,Dl+Ul+Dr+Ur,0));
    G = Ge[2:(nc+1),2:(nc+1)]

    ## Handle boundary conditions
    if(bc[1]=='p') G[1,nc] = Ge[2,1]         ## Exit at left to re-entry at right
    if(bc[2]=='p') G[nc,1] = Ge[nc+1,nc+2]   ## Exit at right to re-entry at left

    if(bc[1]=='r') G[1,1] <- -G[1,2]         ## Reflect
    if(bc[2]=='r') G[nc,nc] <- -G[nc,nc-1]   

    ## Experimental: Continue beyod boundary by neglecting diffusivity at that boundary
    if(bc[1]=='c') G[1,1] <- - (G[1,2] <- Ur[1])
    if(bc[2]=='c') G[nc,nc] <- - (G[nc,nc-1] <- Ul[nc])

    if(bc[1]=='e') G <- Ge;
    
    if(sparse) return(G) else return(as.matrix(G))
}
