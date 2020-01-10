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

#' Compute the stationary distribution for a Continuous Time Markov Chain
#'
#' @description
#' StationaryDistribution compute the stationary distribution phi for a Continuous Time
#' Markov Chain given by a generator G, which is given by the equations
#'
#'   phi %*% G == 0, sum(phi) == 1
#'
#' It is assumed that these two equations specify phi uniquely, i.e. that G is the generator
#' of an ergodic chain.
#'
#' @param G Generator of a CTMC: A quadratic matrix with non-negative off-diagonal elements and zero row sums
#'
#' @examples
#' G <- arrayc(c(-2,1,2,-1),c(2,2))
#' phi <- StationaryDistribution(G)
#'
#' @export
StationaryDistribution <- function(G)
{
    require(Matrix)
    
    ## Pad with ones to set up for solution
    Gpad <- rbind(cbind(G,rep(1,nrow(G))),c(rep(1,ncol(G)),0))

    b <- sparseMatrix(i = ncol(G)+1,j=1,x=1)
    tG <- Matrix::t(Gpad)

    phiE <- as.numeric(Matrix::solve(tG,b))

    phi <- phiE[1:ncol(G)]
    return(phi)
}

#' Convert cell probabilities to (average) probability densities
#'
#' @description
#' prob2pdf takes a vector of probabilities and a spatial grid (given by end-points and interfaces)
#' and computes the average p.d.f. in each grid cell
#'
#' @param phi Vector of probabilities
#' @param xgrid Spatial grid
#' @param ygrid Optional grid in y-direction (see details)
#'
#' @details
#' length(phi) must equal length(xgrid)-1, if ygrid is not given, or
#' (length(xgrid)-1)*(length(ygrid)-1) if ygrid is given
#'
#' xgrid (and ygrid, if given) must be strictly increasing
#'
#' @examples
#' xgrid <- seq(-2,2,length=10)^3
#' phi <- diff(pnorm(xgrid))
#' f <- prob2pdf(phi)
#' plot(head(xgrid,-1),f,type="S")
#'
#' @export
prob2pdf <- function(phi,xgrid=xgrid,ygrid=ygrid)
    return(phi/rep(diff(ygrid),length(xgrid)-1)/rep(diff(xgrid),rep(length(ygrid)-1,length(xgrid)-1)))


#' Convert a tabulated function on the plane to a vector
#'
#' @description
#' pack.field takes an array (e.g., of probabilities) and produces a vector. 
#'
#' @param phi Array of probabilities for a 2D grid
#'
#' @export
pack.field <- function(C) as.numeric(C)


#' Convert a vector to a function on the plane
#'
#' @description
#' unpack.field takes a vector and produces a 2D array. 
#'
#' @param phi Tabulated function
#' @param nx Number of grid cells in the x direction
#' @param ny Number of grid cells in the y direction
#' 
#' @details
#' length(phi) must equal nx * ny. 
#' @export
unpack.field <- function(C,nx=length(xc),ny=length(yc)) array(C,c(ny,nx))


#' Compute generator for a 2D advection-diffusion equation
#' @export
fvade2d <- function(ux,uy,Dx,Dy,xgrid,ygrid,sparse=TRUE)
{
    require(Matrix)
    
    ## Number of grid cells
    ncx = length(xgrid)-1;
    dx = diff(xgrid);
    ncy = length(ygrid)-1;
    dy = diff(ygrid);

    xc = 0.5*(head(xgrid,-1)+tail(xgrid,-1))
    yc = 0.5*(head(ygrid,-1)+tail(ygrid,-1))


    ## Ordering of cells
    ## Start with the (x,y) plane. Cell #1 is in the lower left corner. Cell #2 is just above it (N).
    ## Index i corresponds to the y coordinate, index j to the x
    ## Remember that we extend the grid to include a boundary surrounding it,
    ## so i and j run from 0 to length(yc)+1 and length(xc)+1

    ## Cell number for a cell, as well as for its N,S,E,W neighbours. Note that
    ## neighbours assume that these are interior cells (not boundary cells)
    C <- function(i,j) i+j*(length(yc)+2)+1
    N <- function(i,j) C(i,j)+1
    E <- function(i,j) C(i,j)+length(yc)+2
    W <- function(i,j) C(i,j)-length(yc)-2
    S <- function(i,j) C(i,j)-1

    ind2ij <- function(k)
    {
        j <- (k-1) %/% (length(yc)+2) 
        i <- k - j * (length(yc) +2) - 1
        return(c(i,j))
    }

    ## Test hopping back and forth
    ## for(i in 0:(length(yc)+1))for(j in 0:(length(xc)+1)) { ij <- ind2ij(C(i,j)) ; print(c(i,j,ij-c(i,j))) }


    GDentries <- array(0,c(length(xc)*length(yc)*4,3))
    GAentries <- array(0,c(length(xc)*length(yc)*4,3))

    dx <- diff(xgrid)
    dxc <- diff(xc)
    dxc <- c(dxc[1],dxc,tail(dxc,1))
    
    dy <- diff(ygrid)
    dyc <- diff(yc)
    dyc <- c(dyc[1],dyc,tail(dyc,1))

    k <- 0
    l <- 0
    for(i in 1:length(yc))
        for(j in 1:length(xc))    
        {
            ## Diffusive flux N
            k <- k + 1 
            GDentries[k,1] <- C(i,j)
            GDentries[k,2] <- N(i,j)
            GDentries[k,3] <- Dy(xc[j],ygrid[i+1])/dy[i]/dyc[i+1]

            ## Diffusive flux S
            k <- k + 1 
            GDentries[k,1] <- C(i,j)
            GDentries[k,2] <- S(i,j)
            GDentries[k,3] <- Dy(xc[j],ygrid[i])/dy[i]/dyc[i]

            ## Diffusive flux E
            k <- k + 1 
            GDentries[k,1] <- C(i,j)
            GDentries[k,2] <- E(i,j)
            GDentries[k,3] <- Dx(xgrid[j+1],yc[i])/dx[j]/dxc[j+1]

            ## Diffusive flux W
            k <- k + 1 
            GDentries[k,1] <- C(i,j)
            GDentries[k,2] <- W(i,j)
            GDentries[k,3] <- Dx(xgrid[j],yc[i])/dx[j]/dxc[j]

            ## Advective flux N
            l <- l + 1 
            GAentries[l,1] <- C(i,j)
            GAentries[l,2] <- N(i,j)
            GAentries[l,3] <- max(0,uy(xc[j],ygrid[i+1])/dy[i])

            ## Advective flux S
            l <- l + 1 
            GAentries[l,1] <- C(i,j)
            GAentries[l,2] <- S(i,j)
            GAentries[l,3] <- max(0,-uy(xc[j],ygrid[i])/dy[i])

            ## Advective flux E
            l <- l + 1 
            GAentries[l,1] <- C(i,j)
            GAentries[l,2] <- E(i,j)
            GAentries[l,3] <- max(0,ux(xgrid[j+1],yc[i])/dx[j])

            ## Advective flux W
            l <- l + 1 
            GAentries[l,1] <- C(i,j)
            GAentries[l,2] <- W(i,j)
            GAentries[l,3] <- max(0,-ux(xgrid[j],yc[i])/dx[j])
        }

    Gent <- rbind(GDentries,GAentries)
    
    dims=rep( (length(yc)+2)*(length(xc)+2),2)
    Ge <- sparseMatrix(i=Gent[,1],j=Gent[,2],x=Gent[,3],dims = dims)

    ## Boundary conditions: Reflective, for now
    Iinner <- as.numeric(outer(1:length(yc),1:length(xc),Vectorize(C)))

    G <- Ge[Iinner,Iinner]

    G <- G + Diagonal(n=nrow(G),x=-Matrix::rowSums(G))
}

