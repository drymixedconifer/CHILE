# Core standalone generators; all run with terra and base R.

nlm_random <- function(ncol=256,nrow=256,resolution=1,minval=0,maxval=1,template=NULL,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  r <- .chile_raster(ncol,nrow,resolution,template)
  terra::values(r) <- stats::runif(terra::ncell(r),minval,maxval)
  names(r) <- 'random'; r
}

nlm_planargradient <- function(ncol=256,nrow=256,resolution=1,angle=0,template=NULL) {
  r <- .chile_raster(ncol,nrow,resolution,template)
  xy <- terra::xyFromCell(r,seq_len(terra::ncell(r)))
  z <- xy[,1]*cos(angle*pi/180)+xy[,2]*sin(angle*pi/180)
  z <- (z-min(z))/max(diff(range(z)), .Machine$double.eps)
  terra::values(r) <- z; names(r) <- 'gradient'; r
}

# Exponential covariance, including geometric anisotropy. On the FFT torus,
# the spectral density in 2D is proportional to (1 + k^2)^(-3/2).
nlm_gaussianfield <- function(ncol=256,nrow=256,resolution=1,range=20,sill=1,
                              mean=0,angle=0,range_minor=range,
                              template=NULL,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (range<=0 || range_minor<=0 || sill<0) stop('Invalid Gaussian field parameters')
  r <- .chile_raster(ncol,nrow,resolution,template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  fx <- matrix(rep(.chile_freq(nc)/terra::res(r)[1],each=nr),nr,nc)
  fy <- matrix(rep(.chile_freq(nr)/terra::res(r)[2],times=nc),nr,nc)
  a <- angle*pi/180
  u <- fx*cos(a)+fy*sin(a)
  v <- -fx*sin(a)+fy*cos(a)
  amp <- (1+(2*pi)^2*((range*u)^2+(range_minor*v)^2))^(-3/4)
  .chile_out(mean+sqrt(sill)*.chile_fft_field(nr,nc,amp),r,'gaussian_exponential')
}

# Matern approximation using its spectral density in two dimensions.
nlm_gaussianfield_matern <- function(ncol=256,nrow=256,resolution=1,range=20,sill=1,mean=0,nu=1,angle=0,range_minor=range,template=NULL,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (range<=0 || range_minor<=0 || nu<=0 || sill<0) stop('Invalid Matern parameters')
  r <- .chile_raster(ncol,nrow,resolution,template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  fx <- matrix(rep(.chile_freq(nc)/terra::res(r)[1],each=nr),nr,nc)
  fy <- matrix(rep(.chile_freq(nr)/terra::res(r)[2],times=nc),nr,nc)
  a <- angle*pi/180
  u <- fx*cos(a)+fy*sin(a); v <- -fx*sin(a)+fy*cos(a)
  # Spectrum proportional to (2 nu + (2 pi)^2 scaled frequency squared)^-(nu+1)
  amplitude <- (2*nu+(2*pi)^2*((range*u)^2+(range_minor*v)^2))^(-(nu+1)/2)
  .chile_out(mean+sqrt(sill)*.chile_fft_field(nr,nc,amplitude),r,'matern')
}

nlm_percolation <- function(ncol=256,nrow=256,resolution=1,prob=0.5,template=NULL,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (prob<0 || prob>1) stop('prob must be between zero and one')
  r <- .chile_raster(ncol,nrow,resolution,template)
  terra::values(r) <- stats::rbinom(terra::ncell(r),1,prob)
  names(r) <- 'percolation'; r
}

# Spectral approximation to fractional Brownian surfaces, anchored at one pixel.
nlm_fbm <- function(ncol=256,nrow=256,resolution=1,H=0.7,template=NULL,seed=NULL) {
  if (H<=0 || H>=1) stop('H must be between zero and one')
  r <- nlm_spectral(ncol=ncol,nrow=nrow,resolution=resolution,beta=2*H+2,
                    template=template,seed=seed)
  v <- terra::values(r)[,1]; terra::values(r) <- v-v[1]
  names(r) <- 'fbm_approx'; r
}

# Genuine diamond-square midpoint displacement. Simulates the next enclosing
# 2^k+1 square then crops to the requested raster grid. roughness controls
# the decay of stochastic displacement at successive refinement levels.
nlm_midpoint <- function(size=129,roughness=0.7,ncol=size,nrow=size,
                         resolution=1,template=NULL,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.finite(roughness) || roughness<=0 || roughness>=1)
    stop('roughness must be in (0,1)')
  r <- .chile_raster(ncol,nrow,resolution,template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  side <- 2^ceiling(log2(max(nr,nc)-1))+1L
  side <- max(3L,as.integer(side))
  z <- matrix(NA_real_,side,side)
  z[c(1,side),c(1,side)] <- stats::rnorm(4)
  step <- side-1L
  amplitude <- 1
  while (step>1L) {
    half <- step %/% 2L
    # Diamond step
    for (i in seq.int(1L,side-step,by=step))
      for (j in seq.int(1L,side-step,by=step))
        z[i+half,j+half] <- mean(c(z[i,j],z[i+step,j],
          z[i,j+step],z[i+step,j+step]))+stats::rnorm(1,sd=amplitude)
    # Square step; boundaries are averaged using existing neighbors.
    for (i in seq.int(1L,side,by=half)) {
      jstart <- if (((i-1L) %/% half) %% 2L == 0L) 1L+half else 1L
      for (j in seq.int(jstart,side,by=step)) {
        neighbors <- numeric(0)
        if (i-half>=1L) neighbors <- c(neighbors,z[i-half,j])
        if (i+half<=side) neighbors <- c(neighbors,z[i+half,j])
        if (j-half>=1L) neighbors <- c(neighbors,z[i,j-half])
        if (j+half<=side) neighbors <- c(neighbors,z[i,j+half])
        z[i,j] <- mean(neighbors)+stats::rnorm(1,sd=amplitude)
      }
    }
    amplitude <- amplitude*roughness
    step <- half
  }
  z <- z[seq_len(nr),seq_len(nc),drop=FALSE]
  z <- z-min(z)
  if (max(z)>0) z <- z/max(z)
  .chile_out(z,r,'midpoint')
}

# Native ambient algorithms live in R/ambient.R.

nlm_mosaic <- function(ncol=256,nrow=256,resolution=1,npts=20,
                       angle=0,range_major=1,range_minor=1,template=NULL,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (npts<1 || range_major<=0 || range_minor<=0) stop('Invalid mosaic parameters')
  r <- .chile_raster(ncol,nrow,resolution,template)
  pts <- .chile_xy(r,npts)
  cells <- terra::xyFromCell(r,seq_len(terra::ncell(r)))
  a <- angle*pi/180
  best <- rep(Inf,nrow(cells)); id <- integer(nrow(cells))
  for (i in seq_len(npts)) {
    dx <- cells[,1]-pts[i,1]; dy <- cells[,2]-pts[i,2]
    u <- (dx*cos(a)+dy*sin(a))/range_major
    v <- (-dx*sin(a)+dy*cos(a))/range_minor
    d <- u^2+v^2
    hit <- d<best; best[hit] <- d[hit]; id[hit] <- i
  }
  terra::values(r) <- id; names(r) <- 'mosaic'; r
}

nlm_mosaiccluster <- function(ncol=256,nrow=256,resolution=1,npts=20,cluster_sd=5,template=NULL,seed=NULL) {
  # cluster_sd had no effect in the original; this alias is intentionally explicit.
  nlm_mosaic(ncol,nrow,resolution,npts,template=template,seed=seed)
}
