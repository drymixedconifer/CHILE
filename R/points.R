# Point-process generators. All positions are continuous XY map coordinates.
# n fixes the count when supplied; n=NULL draws process-specific counts.

.chile_xy <- function(r, n) {
  e <- terra::ext(r)
  cbind(x = stats::runif(n, e[1], e[2]),
        y = stats::runif(n, e[3], e[4]))
}

.chile_in <- function(xy, r) {
  if (!nrow(xy)) return(logical(0))
  e <- terra::ext(r)
  xy[, 1] >= e[1] & xy[, 1] < e[2] &
    xy[, 2] >= e[3] & xy[, 2] < e[4]
}

.chile_sample_xy <- function(r, n, environment = NULL) {
  if (n == 0) return(matrix(numeric(), 0, 2, dimnames=list(NULL,c('x','y'))))
  if (is.null(environment)) return(.chile_xy(r, n))
  if (!terra::compareGeom(r, environment, stopOnError = FALSE))
    stop('environment grid does not match output grid')
  v <- terra::values(environment)[, 1]
  v[!is.finite(v)] <- 0
  if (any(v < 0)) stop('environment weights must be nonnegative')
  if (!any(v > 0)) stop('environment has no positive weights')
  cell <- sample.int(terra::ncell(r), n, replace = TRUE, prob = v)
  xy <- terra::xyFromCell(r, cell)
  xy[, 1] <- xy[, 1] + stats::runif(n, -0.5, 0.5) * terra::res(r)[1]
  xy[, 2] <- xy[, 2] + stats::runif(n, -0.5, 0.5) * terra::res(r)[2]
  colnames(xy) <- c('x','y')
  xy
}

.chile_thin <- function(xy, d, max_n = Inf) {
  if (nrow(xy) < 2 || d <= 0) return(xy[seq_len(min(nrow(xy), max_n)), , drop=FALSE])
  chosen <- integer(0)
  for (i in sample.int(nrow(xy))) {
    if (!length(chosen) || all(rowSums((xy[chosen, , drop=FALSE] -
      matrix(xy[i, ], length(chosen), 2, byrow=TRUE))^2) >= d^2))
      chosen <- c(chosen, i)
    if (length(chosen) >= max_n) break
  }
  xy[chosen, , drop=FALSE]
}

#' Simulate spatial point processes
#'
#' Standalone point-process simulator. Returns subcell XY points and cell counts.
#' `n` requests an exact count for most models; use `n = NULL` with `intensity`
#' for an unconditional Poisson count or naturally random cluster offspring counts.
#' `strauss` is an approximate finite-window birth/death sampler, not an exact
#' perfect sampler. `inhibition` is sequential inhibition, not a Gibbs hard-core model.
#'
#' @param method One of `binomial`, `poisson`, `inhomogeneous`, `thomas`,
#'   `matern_cluster`, `inhibition`, `strauss`, `lgcp`, `grid`,
#'   `jittered`, `anisotropic_cluster`, `mixture`.
#' @param n Requested point count; `NULL` permits model-derived counts in
#'   Poisson, cluster and LGCP modes. Some methods condition on `n`.
#' @param intensity Mean point density per map-unit squared for Poisson modes.
#' @param environment Optional nonnegative raster for intensity weighting.
#' @param parents Number of cluster parents.
#' @param offspring Expected offspring per parent if `n=NULL`.
#' @param cluster_sd Thomas cluster standard deviation in map units.
#' @param radius Cluster radius for Matern clusters; interaction radius for Strauss.
#' @param min_dist Sequential-inhibition minimum separation in map units.
#' @param interaction Strauss pair interaction parameter; below 1 is repulsion.
#' @param range_major,range_minor Major/minor scales for anisotropic clusters or LGCP.
#' @param angle Orientation in degrees counterclockwise from x-axis.
#' @param spacing Grid point spacing in map units; derived from `n` if `NULL`.
#' @param jitter Fraction of grid spacing used for random displacement.
#' @param mixture Vector of component methods (for `method="mixture"`).
#' @param weights Relative proportions of mixture components.
#' @param max_attempts Upper limit on rejection/Metropolis proposals.
#' @param template Optional terra raster defining output grid.
#' @param ncol,nrow,resolution Output dimensions when template is absent.
#' @param seed Optional random seed.
#' @return List containing `raster` (point-count SpatRaster), `points` (data.frame
#'   with continuous XY), and `parameters` (requested model settings).
#' @examples
#' p <- nlm_points(method="thomas", n=200, seed=1)
#' plot(p$raster)
#' plot(p$points$x, p$points$y, asp=1)
#' @export
nlm_points <- function(ncol=256, nrow=256, resolution=1,
                       method=c('binomial','poisson','inhomogeneous','thomas',
                         'matern_cluster','inhibition','strauss','lgcp','grid',
                         'jittered','anisotropic_cluster','mixture'),
                       n=100, intensity=NULL, environment=NULL,
                       parents=12, offspring=10, cluster_sd=12,
                       radius=15, min_dist=8, interaction=0.5,
                       range_major=35, range_minor=8, angle=0,
                       spacing=NULL, jitter=0.3,
                       mixture=c('thomas','inhibition'), weights=NULL,
                       max_attempts=100000, template=NULL, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  method <- match.arg(method)
  r <- .chile_raster(ncol, nrow, resolution, template)
  e <- terra::ext(r)
  area <- (e[2]-e[1]) * (e[4]-e[3])
  if (!is.null(n) && (length(n)!=1 || n < 0 || n != floor(n)))
    stop('n must be a nonnegative integer or NULL')
  if (parents < 1 || offspring <= 0 || cluster_sd < 0 || radius <= 0 ||
      min_dist < 0 || range_major <= 0 || range_minor <= 0 ||
      jitter < 0 || jitter > 1 || max_attempts < 1)
    stop('Invalid point-process parameter')
  if (!is.null(intensity) && (!is.finite(intensity) || intensity < 0))
    stop('intensity must be nonnegative')
  if (method == 'strauss' && (interaction < 0 || interaction > 1))
    stop('strauss interaction must be between 0 and 1')
  count <- function() if (!is.null(n)) n else
    stats::rpois(1, if (is.null(intensity)) 100 else intensity * area)
  empty <- matrix(numeric(0), ncol=2, dimnames=list(NULL,c('x','y')))
  xy <- empty

  if (method %in% c('binomial','poisson','inhomogeneous')) {
    if (method == 'poisson' && !is.null(n)) {
      # Count-conditioned Poisson equals binomial point process.
      # To draw Poisson count supply n=NULL.
    }
    xy <- .chile_sample_xy(r, count(), environment)

  } else if (method %in% c('thomas','matern_cluster','anisotropic_cluster')) {
    nparents <- if (is.null(n)) stats::rpois(1, parents) else parents
    pp <- .chile_xy(r, nparents)
    k <- if (is.null(n)) stats::rpois(nparents, offspring) else {
      z <- integer(nparents)
      if (nparents > 0 && n > 0)
        z <- as.vector(stats::rmultinom(1, n, rep(1, nparents)))
      z
    }
    if (sum(k) > 0) {
      parent <- rep(seq_len(nparents), times=k)
      if (method == 'matern_cluster') {
        theta <- stats::runif(length(parent), 0, 2*pi)
        rad <- radius * sqrt(stats::runif(length(parent)))
        offset <- cbind(rad*cos(theta), rad*sin(theta))
      } else {
        offset <- cbind(stats::rnorm(length(parent)),stats::rnorm(length(parent)))
        if (method == 'anisotropic_cluster') {
          a <- angle*pi/180
          xx <- offset[,1]*range_major
          yy <- offset[,2]*range_minor
          offset <- cbind(xx*cos(a)-yy*sin(a),xx*sin(a)+yy*cos(a))
        } else offset <- offset*cluster_sd
      }
      xy <- pp[parent, , drop=FALSE] + offset
      xy <- xy[.chile_in(xy,r),,drop=FALSE]
      colnames(xy) <- c('x','y')
      if (!is.null(n) && nrow(xy) < n) {
        # Edge correction: resample candidate offspring until count reached.
        attempts <- 0L
        while (nrow(xy) < n && attempts < max_attempts) {
          q <- pp[sample.int(nparents,1),]
          if (method == 'matern_cluster') {
            a <- stats::runif(1,0,2*pi); d <- radius*sqrt(stats::runif(1))
            q <- q + c(cos(a),sin(a))*d
          } else {
            uv <- stats::rnorm(2)
            if (method == 'anisotropic_cluster') {
              a <- angle*pi/180
              u <- uv[1]*range_major; v <- uv[2]*range_minor
              q <- q + c(u*cos(a)-v*sin(a),u*sin(a)+v*cos(a))
            } else q <- q + uv*cluster_sd
          }
          if (.chile_in(matrix(q,1,2),r)) xy <- rbind(xy,q)
          attempts <- attempts+1L
        }
        if (nrow(xy) < n) warning('Cluster edge correction reached max_attempts')
      }
    }

  } else if (method == 'inhibition') {
    target <- count()
    xy <- empty
    attempts <- 0L
    while (nrow(xy) < target && attempts < max_attempts) {
      q <- .chile_sample_xy(r,1,environment)
      if (!nrow(xy) || all(rowSums((xy-matrix(q,nrow(xy),2,byrow=TRUE))^2)>=min_dist^2))
        xy <- rbind(xy,q)
      attempts <- attempts+1L
    }
    if (nrow(xy) < target) warning('Inhibition reached max_attempts; returned fewer points')

  } else if (method == 'strauss') {
    target <- count()
    xy <- .chile_sample_xy(r,target,environment)
    # Fixed-N Metropolis relocation conditional on target count.
    for (iter in seq_len(min(max_attempts, 50L*max(1L,target)))) {
      if (!target) break
      i <- sample.int(target,1)
      q <- .chile_sample_xy(r,1,environment)
      other <- xy[-i,,drop=FALSE]
      old_pairs <- if (nrow(other)) sum(rowSums((other-matrix(xy[i,],nrow(other),2,byrow=TRUE))^2)<radius^2) else 0
      new_pairs <- if (nrow(other)) sum(rowSums((other-matrix(q,nrow(other),2,byrow=TRUE))^2)<radius^2) else 0
      accept <- if (interaction==0) new_pairs <= old_pairs else
        stats::runif(1) < min(1, interaction^(new_pairs-old_pairs))
      if (accept) xy[i,] <- q
    }

  } else if (method == 'lgcp') {
    latent <- nlm_anisotropic(range_major=range_major,range_minor=range_minor,
       angle=angle,template=r)
    z <- terra::values(latent)[,1]
    v <- exp(pmin(20, z - max(z)))
    if (!is.null(environment)) {
      if (!terra::compareGeom(r,environment,stopOnError=FALSE)) stop('environment mismatch')
      ev <- terra::values(environment)[,1]
      ev[!is.finite(ev)] <- 0
      if (any(ev<0)) stop('environment must be nonnegative')
      v <- v*ev
    }
    if (!any(v>0)) stop('LGCP intensity is zero everywhere')
    if (is.null(n)) {
      target_mean <- if (is.null(intensity)) 100 else intensity*area
      lambda <- v/sum(v)*target_mean
      cells <- rep(seq_along(v),stats::rpois(length(v),lambda))
      xy <- terra::xyFromCell(r,cells)
      if (length(cells)) {
        xy[,1] <- xy[,1]+stats::runif(length(cells),-0.5,0.5)*terra::res(r)[1]
        xy[,2] <- xy[,2]+stats::runif(length(cells),-0.5,0.5)*terra::res(r)[2]
      }
      colnames(xy) <- c('x','y')
    } else {
      w <- r; terra::values(w) <- v
      xy <- .chile_sample_xy(r,n,w)
    }

  } else if (method %in% c('grid','jittered')) {
    target <- count()
    if (is.null(spacing)) spacing <- sqrt(area/max(target,1))
    if (spacing<=0) stop('spacing must be positive')
    xx <- seq(e[1]+spacing/2,e[2],by=spacing)
    yy <- seq(e[3]+spacing/2,e[4],by=spacing)
    if (length(xx) && length(yy)) {
      xy <- as.matrix(expand.grid(x=xx,y=yy))
      a <- angle*pi/180
      center <- c((e[1]+e[2])/2,(e[3]+e[4])/2)
      uv <- sweep(xy,2,center)
      xy <- cbind(uv[,1]*cos(a)-uv[,2]*sin(a),
                  uv[,1]*sin(a)+uv[,2]*cos(a))
      xy <- sweep(xy,2,center,'+')
      if (method=='jittered') xy <- xy + matrix(stats::runif(length(xy),-jitter*spacing/2,jitter*spacing/2),ncol=2)
      xy <- xy[.chile_in(xy,r),,drop=FALSE]
      if (!is.null(n) && nrow(xy)>n) xy <- xy[sample.int(nrow(xy),n),,drop=FALSE]
    }
    colnames(xy) <- c('x','y')

  } else if (method == 'mixture') {
    if (!length(mixture) || any(mixture=='mixture')) stop('Invalid mixture')
    target <- count()
    if (is.null(weights)) weights <- rep(1,length(mixture))
    if (length(weights)!=length(mixture) || any(weights<0) || !any(weights>0)) stop('Invalid mixture weights')
    ks <- as.vector(stats::rmultinom(1,target,weights))
    chunks <- lapply(seq_along(mixture),function(j) {
      if (ks[j]==0) return(empty)
      z <- nlm_points(method=mixture[j],n=ks[j],template=r,environment=environment,
        parents=parents,offspring=offspring,cluster_sd=cluster_sd,radius=radius,
        min_dist=min_dist,interaction=interaction,range_major=range_major,
        range_minor=range_minor,angle=angle,spacing=spacing,jitter=jitter,
        max_attempts=max_attempts)
      as.matrix(z$points[,c('x','y'),drop=FALSE])
    })
    xy <- do.call(rbind,chunks)
  }

  if (is.null(dim(xy))) xy <- matrix(xy,ncol=2)
  colnames(xy) <- c('x','y')
  xy <- xy[.chile_in(xy,r),,drop=FALSE]
  cells <- if(nrow(xy)) terra::cellFromXY(r,xy) else integer(0)
  counts <- tabulate(cells,nbins=terra::ncell(r))
  terra::values(r) <- counts
  names(r) <- 'point_count'
  pts <- data.frame(id=seq_len(nrow(xy)),x=xy[,1],y=xy[,2])
  list(raster=r,points=pts,parameters=list(method=method,n_requested=n,
    n_realized=nrow(xy),intensity=intensity,angle=angle))
}
