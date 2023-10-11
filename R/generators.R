#' Discretize scalar advection-diffusion equation using finite volumes
#'
#' @name fvade
#'
#' @description
#' fvade discretizes the advection-diffusion equation
#'
#'      dC/dt = -(u*C-D*C')'
#'
#'   on an interval [a,b] using the finite-volume method. A loss term can be included, as in
#'
#'      dC/dt = -(u*C-D*C')' - r*C
#'
#' @param  u function mapping state (numeric scalar) to advective term (numeric scalar)
#' @param  D function mapping state (numeric scalar) to diffusivity (numeric scalar)
#' @param  xgrid The numerical grid. Numeric vector of increasing values, giving cell boundaries
#' @param  bc String indicating boundary conditions. See details.
#' @param  sparse logical indicating if the result should be returned as a sparse matrix
#' @param  diagonals logical indicating if the result should be returned as a list of subdiagonal, diagonal, and superdiagonal
#' @param  r function mapping state (numeric scalar) to mortality/discount rate. Defaults to 0.
#' @return a quadratic matrix, the (sub)generator of the approximating continuous-time Markov chain, with length(xgrid)-1 columns
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
#' phi <- StationaryDistribution(G)         # Find stationary probabilities
#' phi <- phi/dx                            # Convert to densities
#' plot(xc,phi,type="l",xlab="x",ylab="Stationary density")
#'
#' @export
fvade <- function(u,D,xgrid,bc,sparse=TRUE,diagonals=FALSE,r=NULL)
{
    require(Matrix)

    if(length(bc)==1) bc <- rep(bc,2)

    ## Number of grid cells
    nc <- length(xgrid)-1
    dx <- diff(xgrid)

    ## Indeces of left and right cell, with wrap around
    Ileft <- c(nc,1:(nc-1))
    Iright= c(2:nc,1)

    ## Diffusivities at interfaces
    Di <- sapply(xgrid,D)

    ## Diffusion elements of the generator
    Dl <- 2*utils::head(Di,-1)/dx/(dx + dx[Ileft])
    Dr <- 2*utils::tail(Di,-1)/dx/(dx + dx[Iright])

    ## Advection at left and right interfaces
    Ui <- sapply(xgrid,u)
    Uil <- utils::head(Ui,-1)
    Uir <- utils::tail(Ui,-1)

    ## Advection elements of the generator
    Ul <- pmax(-Uil,0)/dx
    Ur <- pmax( Uir,0)/dx

    ## Loss
    if(is.null(r)) Uc <- numeric(nc) else Uc <- sapply(cell.centers(xgrid),r)

    Ge <- bandSparse(nc+2,nc+2,(-1):1,
                     list(c(Dl+Ul,0),
                          -c(0,Dl+Ul+Dr+Ur+Uc,0),
                          c(0,Dr+Ur)))
                     
    G <- Ge[2:(nc+1),2:(nc+1)]

    ## Handle boundary conditions
    if(bc[1]=='p') G[1,nc] <- Ge[2,1]         ## Exit at left to re-entry at right
    if(bc[2]=='p') G[nc,1] <- Ge[nc+1,nc+2]   ## Exit at right to re-entry at left

    if(bc[1]=='r') G[1,1]   <- -Dr[1]  -Ur[1]  -Uc[1]        ## Reflect
    if(bc[2]=='r') G[nc,nc] <- -Dl[nc] -Ul[nc] -Uc[nc]

    ## Experimental: Continue beyond boundary by neglecting diffusivity at that boundary
    if(bc[1]=='c') G[1,1] <- - (G[1,2] <- Ur[1])
    if(bc[2]=='c') G[nc,nc] <- - (G[nc,nc-1] <- Ul[nc])

    if(bc[1]=='e') G <- Ge;

    if(diagonals)
    {
        if(bc[1]=='p') error("With periodic boundary conditions, the system is not tridiagonal")
        return(list(super=G[seq(nc+1,nc*nc-1,nc+1)],diag=G[seq(1,nc*nc,nc+1)],sub=G[seq(2,nc*(nc-1),nc+1)]))
    }

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
#' of an ergodic chain. An error is returned otherwise.
#'
#' @param G Generator of a CTMC: A quadratic matrix with non-negative off-diagonal elements and zero row sums
#' @param tol Tolerance, defaulting to 1e-12, giving the required accuracy on phi %*% G == 0.
#'
#' @return phi, a vector with as many elements as rows in G, containing the stationary distribution.
#'
#' @examples
#' G <- array(c(-2,1,2,-1),c(2,2))
#' phi <- StationaryDistribution(G)
#'
#' @export
StationaryDistribution <- function(G,tol=1e-12)
{
    require(Matrix)

    ## Pad with ones to set up for solution
    Gpad <- rbind(cbind(G,rep(1,nrow(G))),c(rep(1,ncol(G)),0))

    b <- sparseMatrix(i = ncol(G)+1,j=1,x=1)
    tG <- Matrix::t(Gpad)

    phiE <- as.numeric(Matrix::solve(tG,b))

    phi <- phiE[1:ncol(G)]

    lambda <- phiE[ncol(G)+1]
    if(abs(lambda)>tol) stop("The distribution is not stationary. Is G a generator, i.e. singular?")

    return(phi)
}

#' Compute the quasi-stationary distribution for a terminating Continuous Time Markov Chain
#'
#' @description
#' QuasiStationaryDistribution compute the quasi-stationary distribution phi for a Continuous Time
#' Markov Chain given by a sub-generator G, which is given by the equations
#'
#'   phi %*% G == -mu*phi, sum(phi) == 1
#'
#' @param G Sub-generator of a CTMC: A quadratic matrix with non-negative off-diagonal elements and non-positive row sums
#'
#' @return A list containing:
#' vector, a vector with non-negative elements, containing the stationary distribution
#' value, minus the corresponding eigenvalue, i.e., the decay rate.
#'
#' @examples
#' xgrid <- seq(0,1,0.01)
#' G <- fvade(u=function(x)5,D=function(x)1,xgrid=xgrid,bc="a")
#' sol <- QuasiStationaryDistribution(G)
#' plot(cell.centers(xgrid),sol$vector/diff(xgrid),main=paste("Decay rate",sol$value))
#' @export
QuasiStationaryDistribution <- function(G)
{
    evs <- RSpectra::eigs(Matrix::t(G),1,sigma=1e-8)
    pi <- evs$vector
    pi <- Re(pi/sum(pi))

    return(list(value=Re(-evs$value),vector=pi))
}

#' Compute the first transient modes (forward and backward) for a terminating Continuous Time Markov Chain
#'
#' @description
#' TransientModes computes the first transient modes for a Continuous Time
#' Markov Chain generator or sub-generator G. These are given by the equations
#'
#'   phi %*% G == lambda*phi, G %*% psi == lambda*psi
#'
#' By "first" modes, we mean that the eigenvalues have with the largest real parts
#'
#' If G is a generator, then the first mode is the stationary mode. Then the first left mode right mode is
#' trivial, i.e., the constant function, while the first left mode is the stationary distribution. The
#' second value is one divided by the decorrelation time.
#'
#' @param G (Sub-)generator of a CTMC: A quadratic matrix with non-negative off-diagonal elements and non-positive row sums
#' @param k Number of modes. 
#'
#' @return A list containing:
#' phi, a matrix so that each of the k rows contains a transient mode of the forward Kolmogorov equation
#' psi, a matrix so that each of the k columns contains a transient mode of the backward Kolmogorov equation
#' lambda, a vector of the k eigenvalues
#'
#' @examples
#' ## Modes for the Ornstein-Uhlenbeck process
#' xgrid <- seq(-3,3,0.01)
#' xc <- cell.centers(xgrid)
#' G <- fvade(u=function(x)-x,D=function(x)1,xgrid=xgrid,bc="r")
#' modes <- TransientModes(G,3)
#' par(mfrow=c(3,2))
#' for(i in 1:3) {
#'   plot(xc,modes$phi[i,],type="l",main=modes$lambda[i])
#'   plot(xc,modes$psi[,i],type="l",main=modes$lambda[i])
#' }
#' @export
TransientModes <- function(G,k=5)
{
    phi <- t(RSpectra::eigs(Matrix::t(G),k,sigma=1e-8)$vector)
    phi <- phi[k:1,]
    phi[1,] <- phi[1,] / sum(phi[1,])

    evs <- RSpectra::eigs(G,k,sigma=1e-8)
    psi <- evs$vector[,k:1]
    psi[,1] <- psi[,1] / mean(psi[,1])

    lambda <- rev(evs$value)

    return(list(phi=phi,psi=psi,lambda=lambda))
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
#' f <- prob2pdf(phi,xgrid)
#' plot(utils::head(xgrid,-1),f,type="S")
#'
#' @export
prob2pdf <- function(phi,xgrid,ygrid=c(0,1))
    return(phi/rep(diff(ygrid),length(xgrid)-1)/rep(diff(xgrid),rep(length(ygrid)-1,length(xgrid)-1)))


#' Convert a tabulated function on the plane to a vector
#'
#' @description
#' pack.field takes an array (e.g., of probabilities) and produces a vector.
#'
#' @param phi Array of probabilities for a 2D grid
#'
#' @export
pack.field <- function(phi) as.numeric(phi)


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
unpack.field <- function(phi,nx,ny) array(phi,c(ny,nx))

#' Get center of grid cells for a retangular grid
#'
#' @description
#' cell.centers takes a rectangular grid and returns coordinates of the cell centers as a field
#'
#' @param xgrid Interfaces in x-direction
#' @param ygrid Interfaces in y-direction (if 2D)
#'
#' @export
cell.centers <- function(xgrid,ygrid=NULL)
{
    xc <- 0.5*(head(xgrid,-1)+tail(xgrid,-1))

    if(is.null(ygrid)) return(xc)

    yc <- 0.5*(head(ygrid,-1)+tail(ygrid,-1))

    xx <- rep(xc,rep(length(yc),length(xc)))
    yy <- rep(yc,length(xc))

    return(list(x=unpack.field(xx,length(xc),length(yc)),
                y=unpack.field(yy,length(xc),length(yc))))

}

#' Compute generator for a 2D advection-diffusion equation
#'
#' @description
#' fvade2d discretizes the advection-diffusion equation
#'
#'      dC/dt = -div ( u C - D grad C)
#'
#' on a rectangular domain using the finite-volume method. Here, u=(ux,uy) and D=diag(Dx,Dy)
#' @param  ux function mapping state (x,y) to advective term (numeric scalar)
#' @param  uy function mapping state (x,y) to advective term (numeric scalar)
#' @param  Dx function mapping state (x,y) to diffusivity (numeric scalar)
#' @param  Dy function mapping state (x,y) to diffusivity (numeric scalar)
#' @param  xgrid The numerical grid. Numeric vector of increasing values, giving cell boundaries
#' @param  ygrid The numerical grid. Numeric vector of increasing values, giving cell boundaries
#' @param  bc Specification of boundary conditions. See details.
#' @return a quadratic matrix, the generator of the approximating continuous-time Markov chain, with (length(xgrid)-1)*(length(ygrid)-1) columns
#'
#' @details Boundary conditions: bc is a list with elements N,E,S,W. Each element is either
#'   "r": Reflective boundary
#'   "a": Absorbing boundary: Assume an absorbing boundary cell, which is not included
#'   "e": Extend to include an absorbing boundary cell
#'   "p": Periodic. When hitting this boundary, the state is immediately transferred to the opposite boundary, e.g. N->S.
#'
#'  Return value: The function fvade returns a generator (or sub-generator) G of a continuous-time Markov Chain. This chain jumps
#'  between cells defined by xgrid and ygrid. When using the generator to solve the Kolmogorov equations, note that G operates on
#'  probabilities of each cell, not on the probability density in each cell. The distinction is particularly important when the
#'  grid is non-uniform.
#'
#' @examples
#' # Generator of a predator-prey model
#' xi <- seq(0,1.5,0.02)
#' yi <- seq(0,1.6,0.02)
#' xc <- 0.5*(utils::head(xi,-1)+utils::tail(xi,-1))
#' yc <- 0.5*(utils::head(yi,-1)+utils::tail(yi,-1))
#'
#' ux <- function(x,y) x*(1-x)-y*x/(1+x)
#' uy <- function(x,y) y*x/(1+x)-y/3
#' D <- function(x,y) 0.01
#'
#' G <- fvade2d(ux,uy,Dx=D,Dy=D,xi,yi)
#'
#' phiv <- StationaryDistribution(G)
#' phim <- unpack.field(phiv,length(xc),length(yc))
#' image(xi,yi,t(phim))
#'
#' @export
fvade2d <- function(ux,uy,Dx,Dy,xgrid,ygrid,bc=list(N="r",E="r",S="r",W="r"),Dxy=function(x,y)0)
{
    require(Matrix)

    ## Number of grid cells
    ncx <- length(xgrid)-1;
    dx <- diff(xgrid);
    ncy <- length(ygrid)-1;
    dy <- diff(ygrid);

    xc <- 0.5*(utils::head(xgrid,-1)+utils::tail(xgrid,-1))
    yc <- 0.5*(utils::head(ygrid,-1)+utils::tail(ygrid,-1))

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

    GDentries <- array(0,c(length(xc)*length(yc)*4,3))
    GAentries <- array(0,c(length(xc)*length(yc)*4,3))

    dx <- diff(xgrid)
    dxc <- diff(xc)
    dxc <- c(dxc[1],dxc,utils::tail(dxc,1))

    dy <- diff(ygrid)
    dyc <- diff(yc)
    dyc <- c(dyc[1],dyc,utils::tail(dyc,1))

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
    Ge <- Ge + Diagonal(n=nrow(Ge),x=-Matrix::rowSums(Ge))

    ## Handle N,S,E,W boundaries.
    ## First S (given by i==1)
    SC <- C(1,1:ncx)
    SS <- S(1,1:ncx)

    NC <- C(ncy,1:ncx)
    NN <- N(ncy,1:ncx)

    WC <- C(1:ncy,1)
    WW <- W(1:ncy,1)

    EC <- C(1:ncy,ncx)
    EE <- E(1:ncy,ncx)

    ## Which cells to include in the generator that is returned?
    ## The most common choice is only interior cells
    ilo <- 1
    ihi <- ncy

    jlo <- 1
    jhi <- ncx

    ## Reimplementation of extracting diagonal due to weird error in diag
    my.diag <- function(A) as.numeric(A[seq(1,(n<-nrow(A))^2,n+1)])

    ## Process South boundary
    if(bc$S == 'r')
    {
        Ge[SC,SC] <- Ge[SC,SC] + Diagonal(n=ncx,x=my.diag(Ge[SC,SS]))
    }

    if(bc$S == 'e')  ilo <- 0

    if(bc$S == 'p')
    {
        Ge[SC,NC] <- Ge[SC,NC] + Diagonal(n=ncx,x=my.diag(Ge[SC,SS]))
    }

    ## Process North boundary
    if(bc$N == 'r')
    {
        Ge[NC,NC] <- Ge[NC,NC] + Diagonal(n=ncx,x=my.diag(Ge[NC,NN]))
    }

    if(bc$N == 'e')  ihi <- ncy+1

    if(bc$N == 'p')
    {
        Ge[NC,SC] <- Ge[NC,SC] + Diagonal(n=ncx,x=my.diag(Ge[NC,NN]))
    }

    ## Process West boundary
    if(bc$W == 'r')
    {
        Ge[WC,WC] <- Ge[WC,WC] + Diagonal(n=ncy,x=my.diag(Ge[WC,WW]))
    }

    if(bc$W == 'e')  jlo <- 0

    if(bc$W == 'p')
    {
        Ge[WC,EC] <- Ge[WC,EC] + Diagonal(n=ncy,x=my.diag(Ge[WC,WW]))
    }

    ## Process East boundary
    if(bc$E == 'r')
    {
        Ge[EC,EC] <- Ge[EC,EC] + Diagonal(n=ncy,x=my.diag(Ge[EC,EE]))
    }

    if(bc$E == 'e')  jhi <- ncx + 1

    if(bc$E == 'p')
    {
        Ge[EC,WC] <- Ge[EC,WC] + Diagonal(n=ncy,x=my.diag(Ge[EC,EE]))
    }

    ## Extract those cells that should be included
    Iinner <- as.numeric(outer(ilo:ihi,jlo:jhi,Vectorize(C)))

    G <- Ge[Iinner,Iinner]
}
