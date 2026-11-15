# CHILE procedural noise, backed by Thomas Lin Pedersen's ambient package.
# Coordinate transform is shared by noise and patterns; units are raster cells.

.chile_ambient_seed <- function(seed) {
  if (!is.null(seed)) set.seed(seed)
  sample.int(.Machine$integer.max, 1L)
}

.chile_ambient_raw <- function(type, x, y, frequency, seed, interpolator,
                               distance, value, distance_ind, jitter) {
  switch(type,
    perlin = ambient::gen_perlin(x,y,frequency=frequency,seed=seed,interpolator=interpolator),
    simplex = ambient::gen_simplex(x,y,frequency=frequency,seed=seed),
    cubic = ambient::gen_cubic(x,y,frequency=frequency,seed=seed),
    value = ambient::gen_value(x,y,frequency=frequency,seed=seed,interpolator=interpolator),
    white = ambient::gen_white(x,y,frequency=frequency,seed=seed),
    worley = ambient::gen_worley(x,y,frequency=frequency,seed=seed,
       distance=distance,value=value,distance_ind=distance_ind,jitter=jitter),
    waves = ambient::gen_waves(x,y,frequency=frequency),
    spheres = ambient::gen_spheres(x,y,frequency=frequency),
    checkerboard = ambient::gen_checkerboard(x,y,frequency=frequency)
  )
}

#' Native ambient procedural fields
#'
#' A single-raster wrapper around ambient's C++ noise generators. No external
#' environment is required. All frequencies and warp amplitudes are measured
#' in raster-cell coordinates; angle is degrees counterclockwise from x.
#' Supply template to match an existing raster grid. A seed reproduces the field.
#'
#' @param ncol,nrow Grid size without template.
#' @param resolution Cell size without template.
#' @param extent,crs Optional extent and coordinate reference system.
#' @param type One of perlin, simplex, cubic, value, white, worley, blue,
#'   waves, spheres, checkerboard.
#' @param frequency Lattice frequency in inverse cells (e.g., 0.02 means ~50 cells).
#' @param fractal NULL/none, fbm, billow or rigid-multi (rigidmulti accepted).
#' @param octaves,lacunarity,gain Fractal layers, frequency multiplier, amplitude multiplier.
#' @param perturbation none, normal, or fractal coordinate warping.
#' @param perturbation_amplitude Warp displacement in cells (before anisotropic transform).
#' @param angle Major-axis orientation (degrees CCW from east).
#' @param anisotropy Major-to-minor scale ratio; at least 1.
#' @param interpolator linear, hermite, or quintic for Perlin/value noise.
#' @param distance Worley distance: euclidean, manhattan, natural.
#' @param value Worley result: cell, distance, distance2, distance2add,
#'   distance2sub, distance2mul, distance2div.
#' @param distance_ind Two nearest-feature indices for Worley.
#' @param jitter Worley feature displacement.
#' @param blue_sd,blue_seed_frac Blue noise parameters; expensive on large grids.
#' @param normalize Whether to rescale to rescale_range.
#' @param rescale_range Min and max for optional normalization.
#' @param template Optional terra raster template.
#' @param seed Optional integer random seed.
#' @return Single-layer terra SpatRaster.
#' @examples
#' x <- nlm_ambient(type="worley",value="distance2sub",seed=4)
#' y <- nlm_ambient(type="simplex",fractal="billow",angle=35,
#'                  anisotropy=5,perturbation="fractal",seed=4)
#' @export
nlm_ambient <- function(ncol=256,nrow=256,resolution=1,extent=NULL,crs='',
                        type='simplex',frequency=0.02,fractal=NULL,octaves=4,
                        lacunarity=2,gain=0.5,perturbation='none',
                        perturbation_amplitude=1,angle=0,anisotropy=1,
                        interpolator='quintic',distance='euclidean',value='cell',
                        distance_ind=c(1,2),jitter=0.45,
                        blue_sd=10,blue_seed_frac=0.1,
                        normalize=TRUE,rescale_range=c(0,1),
                        template=NULL,seed=NULL) {
  if (!requireNamespace('ambient',quietly=TRUE))
    stop('Install ambient: install.packages("ambient")')
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(extent) && is.null(template))
    template <- terra::rast(nrows=nrow,ncols=ncol,ext=extent,crs=crs)
  r <- .chile_raster(ncol,nrow,resolution,template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  types <- c('perlin','simplex','cubic','value','white','worley','blue',
             'waves','spheres','checkerboard')
  if (length(type)!=1L || !type %in% types)
    stop('type must be one of: ',paste(types,collapse=', '))
  if (identical(fractal,'rigidmulti')) fractal <- 'rigid-multi'
  if (is.null(fractal)) fractal <- 'none'
  if (!fractal %in% c('none','fbm','billow','rigid-multi')) stop('Unknown fractal')
  if (!perturbation %in% c('none','normal','fractal')) stop('Unknown perturbation')
  if (length(frequency)!=1L || !is.finite(frequency) || frequency<=0 ||
      length(anisotropy)!=1L || !is.finite(anisotropy) || anisotropy<1 ||
      !is.finite(angle) || !is.finite(perturbation_amplitude) ||
      perturbation_amplitude<0 || !is.finite(lacunarity) || lacunarity<=0 ||
      !is.finite(gain) || gain<0 || !is.finite(octaves) ||
      octaves<1 || octaves!=as.integer(octaves)) stop('Invalid noise parameters')
  if (!interpolator %in% c('linear','hermite','quintic')) stop('Unknown interpolator')
  if (!distance %in% c('euclidean','manhattan','natural')) stop('Unknown Worley distance')
  if (!value %in% c('cell','distance','distance2','distance2add',
                     'distance2sub','distance2mul','distance2div'))
    stop('Unknown Worley value')
  if (length(distance_ind)!=2L || any(!is.finite(distance_ind)))
    stop('distance_ind needs two finite indices')
  if (length(rescale_range)!=2L || any(!is.finite(rescale_range)) ||
      rescale_range[1]>=rescale_range[2]) stop('Invalid rescale_range')

  if (type=='blue') {
    if (fractal!='none' || perturbation!='none' || anisotropy!=1 || angle!=0)
      stop('blue noise does not support fractal, warp, or anisotropy')
    if (nr*nc>65536L) stop('blue noise limited to 65,536 cells (slow algorithm)')
    z <- ambient::noise_blue(dim=c(nr,nc),sd=blue_sd,seed_frac=blue_seed_frac)
    v <- as.vector(t(z))
  } else {
    # Coordinates in cells; x increases right and y increases upward.
    x <- rep(seq_len(nc)-(nc+1)/2,times=nr)
    y <- rep(((nr+1)/2-seq_len(nr)),each=nc)

    # Warp in original grid coordinates, so amplitude is independent of angle.
    if (perturbation!='none' && perturbation_amplitude>0) {
      nwarp <- if (perturbation=='fractal') 3L else 1L
      wx <- wy <- numeric(length(x))
      for (k in seq_len(nwarp)) {
        freq <- 0.015*2^(k-1L)
        weight <- 0.5^(k-1L)
        wx <- wx + weight*ambient::gen_simplex(x,y,frequency=freq,
                 seed=.chile_ambient_seed(NULL))
        wy <- wy + weight*ambient::gen_simplex(x,y,frequency=freq,
                 seed=.chile_ambient_seed(NULL))
      }
      x <- x+perturbation_amplitude*wx
      y <- y+perturbation_amplitude*wy
    }
    a <- angle*pi/180
    u <- (x*cos(a)+y*sin(a))/anisotropy
    w <- -x*sin(a)+y*cos(a)

    no_fractal <- c('white','waves','spheres','checkerboard')
    if (type %in% no_fractal && fractal!='none')
      stop('Fractal octaves are unsupported for ',type)
    no <- if (fractal=='none') 1L else as.integer(octaves)
    v <- numeric(length(x))
    for (k in seq_len(no)) {
      z <- .chile_ambient_raw(type,u,w,frequency*lacunarity^(k-1L),
             .chile_ambient_seed(NULL),interpolator,distance,value,distance_ind,jitter)
      if (fractal=='billow') z <- 2*abs(z)-1
      if (fractal=='rigid-multi') z <- 1-abs(z)
      v <- v+gain^(k-1L)*z
    }
  }
  if (normalize) {
    limits <- range(v,finite=TRUE)
    if (diff(limits)>0) v <- rescale_range[1]+(v-limits[1]) /
                                     diff(limits)*diff(rescale_range)
    else v[] <- mean(rescale_range)
  }
  terra::values(r) <- v
  names(r) <- paste0('ambient_',type)
  r
}

#' Vector field derived from an ambient surface
#'
#' A gradient points uphill. A 2D curl is its 90-degree rotation,
#' producing a numerically divergence-free interior field under the
#' centered-difference stencil. Edge derivatives are one-sided.
#' @param mode gradient or curl.
#' @param ... Arguments passed to nlm_ambient().
#' @return Two-layer SpatRaster named dx and dy; vectors are in map units.
#' @export
nlm_ambient_flow <- function(mode='curl',...) {
  mode <- match.arg(mode,c('gradient','curl'))
  r <- nlm_ambient(...)
  z <- matrix(terra::values(r)[,1],nrow(r),ncol(r),byrow=TRUE)
  nr <- nrow(z); nc <- ncol(z)
  if (nr<2 || nc<2) stop('Flow needs at least two rows and columns')
  rx <- terra::res(r)[1]; ry <- terra::res(r)[2]
  dx <- cbind((z[,2]-z[,1])/rx,
              if(nc>2) (z[,3:nc,drop=FALSE]-z[,1:(nc-2),drop=FALSE])/(2*rx),
              (z[,nc]-z[,nc-1])/rx)
  dy <- rbind((z[1,]-z[2,])/ry,
              if(nr>2) (z[1:(nr-2),,drop=FALSE]-z[3:nr,,drop=FALSE])/(2*ry),
              (z[nr-1,]-z[nr,])/ry)
  if (mode=='curl') {
    tmp <- dx
    dx <- dy
    dy <- -tmp
  }
  out <- c(r,r)
  terra::values(out) <- cbind(as.vector(t(dx)),as.vector(t(dy)))
  names(out) <- c('dx','dy')
  out
}

#' Compare ambient parameter settings
#'
#' Specify up to two varying arguments to get a small parameter sweep.
#' Settings share a seed but different noise types should not be interpreted
#' as paired stochastic realizations.
#' @param type,fractal,frequency,perturbation_amplitude Parameter candidates.
#' @param ncol,nrow Dimensions of each preview.
#' @param seed Reproducible seed for each panel.
#' @param ... Other scalar nlm_ambient parameters.
#' @return Invisibly returns a named list of rasters.
#' @export
nlm_ambient_show <- function(type='simplex',fractal='none',
                             frequency=0.02,perturbation_amplitude=0,
                             ncol=100,nrow=100,seed=42,...) {
  opts <- list(type=type,fractal=fractal,frequency=frequency,
               perturbation_amplitude=perturbation_amplitude)
  varied <- names(opts)[lengths(opts)>1L]
  if (length(varied)>2L) stop('Vary at most two parameters')
  combos <- expand.grid(opts,stringsAsFactors=FALSE)
  if (nrow(combos)>24L) stop('Maximum 24 panels')
  extra <- list(...)
  if (any(names(extra) %in% names(opts))) stop('Repeated gallery parameter')
  rs <- vector('list',nrow(combos))
  labels <- character(nrow(combos))
  for (i in seq_len(nrow(combos))) {
    p <- as.list(combos[i,,drop=FALSE])
    if (p$perturbation_amplitude>0 && is.null(extra$perturbation))
      p$perturbation <- 'fractal'
    args <- c(list(ncol=ncol,nrow=nrow,seed=seed),p,extra)
    rs[[i]] <- do.call(nlm_ambient,args)
    labels[i] <- if (length(varied))
      paste(paste0(varied,'=',unlist(p[varied])),collapse='  ') else
      paste0('ambient: ',p$type)
  }
  names(rs) <- labels
  old <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(old),add=TRUE)
  graphics::par(mfrow=c(ceiling(nrow(combos)/ceiling(sqrt(nrow(combos)))),
                        ceiling(sqrt(nrow(combos)))),mar=c(2,2,3,1))
  for (i in seq_along(rs)) terra::plot(rs[[i]],main=labels[i],legend=FALSE)
  invisible(rs)
}
