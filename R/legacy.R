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

nlm_gaussianfield <- function(ncol=256,nrow=256,resolution=1,range=20,sill=1,mean=0,angle=0,range_minor=range,template=NULL,seed=NULL) {
  nlm_anisotropic(ncol=ncol,nrow=nrow,resolution=resolution,range_major=range,
                  range_minor=range_minor,angle=angle,mean=mean,sd=sqrt(sill),template=template,seed=seed)
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

nlm_midpoint <- function(size=129,roughness=0.7,ncol=size,nrow=size,resolution=1,template=NULL,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (roughness<=0 || roughness>=1) stop('roughness must be in (0,1)')
  # Multiscale spectral approximation allows arbitrary raster dimensions.
  r <- nlm_multiscale(ncol=ncol,nrow=nrow,resolution=resolution,
    ranges=c(3,10,30,90),weights=roughness^(3:0),template=template)
  names(r) <- 'midpoint_like'; r
}

# Lightweight internal value-noise implementation using multiscale Gaussian fields.
# 'worley' is a distance-to-seed field, and 'white' is independent noise.
nlm_ambient <- function(ncol=256,nrow=256,resolution=1,extent=NULL,crs='',
                        type='perlin',frequency=0.02,fractal=NULL,octaves=4,
                        lacunarity=2,gain=0.5,normalize=TRUE,rescale_range=c(0,1),
                        template=NULL,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(extent) && is.null(template)) {
    template <- terra::rast(nrows=nrow,ncols=ncol,ext=extent,crs=crs)
  }
  r <- .chile_raster(ncol,nrow,resolution,template)
  if (frequency<=0 || octaves<1 || lacunarity<=0 || gain<0) stop('Invalid noise parameters')
  if (type=='white') {
    z <- nlm_random(template=r)
  } else if (type=='worley') {
    nfeat <- max(1,round(terra::ncell(r)*frequency^2))
    z <- nlm_distance(n_features=nfeat,template=r)
  } else if (type %in% c('perlin','simplex','value')) {
    ranges <- 1/(frequency*lacunarity^(seq_len(octaves)-1))
    weights <- gain^(seq_len(octaves)-1)
    z <- nlm_multiscale(ranges=ranges,weights=weights,template=r)
    if (identical(fractal,'billow')) z <- abs(z)
    if (identical(fractal,'rigidmulti')) z <- -abs(z)
  } else stop('Unsupported ambient type')
  v <- terra::values(z)[,1]
  if (normalize) {
    d <- diff(range(v)); if (d>0) v <- (v-min(v))/d else v[] <- mean(rescale_range)
    v <- rescale_range[1]+v*diff(rescale_range)
  }
  terra::values(r) <- v; names(r) <- paste0('ambient_',type); r
}

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
