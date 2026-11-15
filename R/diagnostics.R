# CHILE diagnostics: exploratory, scale-dependent statistics.
# Rasters use planar map distances. Points retain subcell XY locations.

.chile_diag_raster <- function(x) {
  if (is.list(x) && inherits(x$raster, 'SpatRaster')) x <- x$raster
  if (!inherits(x, 'SpatRaster')) stop('Expected a SpatRaster or CHILE result with $raster')
  if (terra::nlyr(x) != 1L) stop('Supply a single raster layer')
  if (terra::is.lonlat(x)) stop('Use a projected raster: distances must be in linear map units')
  x
}

.chile_diag_points <- function(x) {
  if (is.list(x) && !is.null(x$points)) x <- x$points
  if (is.data.frame(x)) {
    if (!all(c('x','y') %in% names(x))) stop('Point table needs x and y columns')
    x <- as.matrix(x[,c('x','y')])
  }
  if (!is.matrix(x) || ncol(x) != 2L || !is.numeric(x))
    stop('Supply CHILE $points, a data.frame with x/y, or two-column XY matrix')
  if (!all(is.finite(x))) stop('Point coordinates must be finite')
  x
}

#' Empirical (directional) semivariogram
#'
#' Randomly samples raster-cell pairs. Distances are in projected map units;
#' direction angles are counterclockwise from +x (modulo 180 degrees).
#' The field is assumed intrinsically stationary for semivariogram interpretation.
#' Returned bin counts give the number of sampled pairs per lag and angle.
#'
#' @param x Continuous, single-layer projected SpatRaster.
#' @param lag Bin width in map units; default approximately two cell widths.
#' @param max_dist Largest lag; defaults to one-quarter of the smaller extent.
#' @param angle Directions in degrees, or NULL for omnidirectional only.
#' @param tolerance Half-width of directional tolerance in degrees.
#' @param pairs Number of randomly sampled cell pairs (before filtering).
#' @param seed Optional seed for reproducibility.
#' @return Data frame: angle, distance, gamma, n_pairs.
#' @examples
#' x <- nlm_anisotropic(ncol=100,nrow=100,seed=1)
#' v <- nlm_variogram(x,angle=c(0,45,90,135),seed=2)
#' plot(gamma~distance,v[v$angle=='0',],type='b')
#' @export
nlm_variogram <- function(x, lag=NULL, max_dist=NULL, angle=NULL,
                          tolerance=22.5, pairs=50000, seed=NULL) {
  r <- .chile_diag_raster(x)
  if (!is.null(seed)) set.seed(seed)
  ext <- terra::ext(r); res <- terra::res(r)
  if (is.null(lag)) lag <- 2*max(res)
  if (is.null(max_dist)) max_dist <- min(ext[2]-ext[1],ext[4]-ext[3])/4
  if (!is.finite(lag) || lag<=0 || !is.finite(max_dist) || max_dist<=lag)
    stop('Require 0 < lag < max_dist')
  if (!is.finite(tolerance) || tolerance<0 || tolerance>90) stop('tolerance must be 0 to 90 degrees')
  if (pairs<1 || pairs!=floor(pairs)) stop('pairs must be a positive integer')
  z <- terra::values(r)[,1]
  cells <- which(is.finite(z))
  if (length(cells)<2) stop('Need two or more valid raster cells')
  i <- sample(cells,pairs,TRUE)
  j <- sample(cells,pairs,TRUE)
  xy1 <- terra::xyFromCell(r,i); xy2 <- terra::xyFromCell(r,j)
  dx <- xy2[,1]-xy1[,1]; dy <- xy2[,2]-xy1[,2]
  d <- sqrt(dx^2+dy^2)
  keep <- d>0 & d<=max_dist
  d <- d[keep]
  theta <- (atan2(dy[keep],dx[keep])*180/pi) %% 180
  sq <- 0.5*(z[i[keep]]-z[j[keep]])^2
  breaks <- seq(0,max_dist,by=lag)
  if (tail(breaks,1) < max_dist) breaks <- c(breaks,max_dist)
  center <- (head(breaks,-1)+tail(breaks,-1))/2
  bins <- findInterval(d,breaks,rightmost.closed=TRUE)
  validbin <- bins>=1 & bins<=length(center)
  dirs <- if (is.null(angle)) NA_real_ else as.numeric(angle)
  if (!length(dirs) || any(!is.finite(dirs) & !is.na(dirs))) stop('Invalid angles')
  out <- lapply(dirs,function(a) {
    selected <- validbin
    if (!is.na(a)) {
      delta <- abs(theta-(a%%180)); delta <- pmin(delta,180-delta)
      selected <- selected & delta<=tolerance
    }
    k <- bins[selected]
    count <- tabulate(k,nbins=length(center))
    total <- numeric(length(center))
    if (length(k)) {
      sums <- rowsum(sq[selected],group=k,reorder=FALSE)
      total[as.integer(rownames(sums))] <- sums[,1]
    }
    data.frame(angle=if (is.na(a)) 'omni' else as.character(a),
               distance=center,gamma=ifelse(count>0,total/pmax(count,1),NA_real_),
               n_pairs=count)
  })
  do.call(rbind,out)
}

#' Ripley's K, transformed L and annular pair correlation
#'
#' Homogeneous finite-window estimator with translation edge correction for a
#' rectangular observation window. K is cumulative; g is annular and may be
#' noisy for narrow radii. L = sqrt(K/pi)-r. Coordinates must use linear units.
#' For nonstationary intensity, this homogeneous statistic can confound
#' variation in intensity with interactions. Limit ~3000 points per pattern.
#'
#' @param x CHILE point result, XY data.frame, or 2-column matrix.
#' @param r Distances in coordinate units. Defaults to 30 evenly spaced distances.
#' @param window Four numbers xmin,xmax,ymin,ymax, or a projected SpatRaster.
#'   CHILE point results automatically use their raster extent.
#' @return Data frame with r, K, L (centered), g, and pair counts.
#' @examples
#' p <- nlm_points(method='thomas',n=150,seed=1)
#' k <- nlm_k(p)
#' plot(L~r,k,type='l')
#' @export
nlm_k <- function(x, r=NULL, window=NULL) {
  if (is.null(window) && is.list(x) && inherits(x$raster,'SpatRaster'))
    window <- x$raster
  xy <- .chile_diag_points(x)
  n <- nrow(xy)
  if (n<2) stop('Need at least two points')
  if (n>3000) stop('For this implementation use <=3000 points or subsample explicitly')
  if (inherits(window,'SpatRaster')) {
    if (terra::is.lonlat(window)) stop('Project coordinates and raster before distance analysis')
    e <- terra::ext(window); window <- c(e[1],e[2],e[3],e[4])
  }
  if (is.null(window)) stop('Provide window=c(xmin,xmax,ymin,ymax) for observed points')
  if (length(window)!=4 || any(!is.finite(window)) || window[2]<=window[1] || window[4]<=window[3])
    stop('Invalid rectangular window')
  if (any(xy[,1]<window[1] | xy[,1]>window[2] | xy[,2]<window[3] | xy[,2]>window[4]))
    stop('Points outside window')
  wx <- window[2]-window[1]; wy <- window[4]-window[3]; A <- wx*wy
  if (is.null(r)) r <- seq(min(wx,wy)/150,min(wx,wy)/4,length.out=30)
  if (!length(r) || any(!is.finite(r)) || any(r<=0) || is.unsorted(r,strictly=TRUE))
    stop('r must be strictly increasing positive distances')
  if (max(r)>=min(wx,wy)) stop('Maximum r must be smaller than both window dimensions')
  d <- stats::dist(xy)
  delta <- as.matrix(d)
  ij <- which(upper.tri(delta),arr.ind=TRUE)
  dx <- abs(xy[ij[,1],1]-xy[ij[,2],1]); dy <- abs(xy[ij[,1],2]-xy[ij[,2],2])
  dp <- delta[ij]
  # Translation overlap area: unbiased for homogeneous stationary processes.
  overlap <- (wx-dx)*(wy-dy)
  weight <- A/overlap
  order_idx <- order(dp)
  dp <- dp[order_idx]; weight <- weight[order_idx]
  cumul <- c(0,cumsum(weight))
  upto <- findInterval(r,dp)
  sumw <- cumul[upto+1L]
  K <- (2*A/(n*(n-1))) * sumw
  L <- sqrt(K/pi)-r
  shell <- c(0,diff(sumw))
  lower <- c(0,head(r,-1))
  shell_area <- pi*(r^2-lower^2)
  g <- (2*A/(n*(n-1)))*shell/shell_area
  data.frame(r=r,K=K,L=L,g=g,n_pairs=upto)
}

#' Connected patch metrics from categorical or thresholded rasters
#'
#' Patches are defined by connected cells of equal class; default 8-neighbor
#' connectivity. The background is NA by default and excluded from metrics.
#' Uses projected linear units; cell sides contribute to perimeter (4 edges).
#' @param x A categorical SpatRaster or CHILE geometry output.
#' @param threshold If supplied, turn a continuous field into binary patches,
#'   selecting values >= threshold; other cells are background NA.
#' @param directions 4 or 8 adjacency for connected patch identities.
#' @return List with per-patch table and class/landscape summary table.
#' @examples
#' b <- nlm_blob(ncol=100,nrow=100,n=15,seed=1)
#' nlm_metrics(b)$summary
#' @export
nlm_metrics <- function(x,threshold=NULL,directions=8) {
  r <- .chile_diag_raster(x)
  if (!directions %in% c(4,8)) stop('directions must be 4 or 8')
  if (!is.null(threshold)) r <- terra::ifel(r>=threshold,1,NA)
  raw <- terra::values(r)[,1]
  vals <- sort(unique(raw[!is.na(raw)]))
  empty_p <- data.frame(class=numeric(),patch=integer(),cells=integer(),area=numeric(),perimeter=numeric(),shape_index=numeric())
  empty_s <- data.frame(class=numeric(),patches=integer(),area=numeric(),mean_area=numeric(),largest_patch_index=numeric(),edge_density=numeric())
  if (!length(vals)) return(list(patches=empty_p,summary=empty_s))
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  rs <- terra::res(r); ca <- prod(rs)
  allp <- vector('list',length(vals)); alls <- vector('list',length(vals))
  for (k in seq_along(vals)) {
    mask <- terra::ifel(r==vals[k],1,NA)
    idr <- terra::patches(mask,directions=directions,zeroAsNA=TRUE)
    z <- matrix(terra::values(idr)[,1],nr,nc,byrow=TRUE)
    pos <- which(is.finite(z))
    if (!length(pos)) next
    ids <- as.integer(z[pos]); n <- max(ids)
    cnt <- tabulate(ids,nbins=n)
    perimeter <- numeric(n)
    # Boundary cells count exposed *sides* (not diagonal neighbor contacts).
    # Edge on top/bottom has width dx; edge left/right has width dy.
    pad <- matrix(0,nr+2L,nc+2L);pad[2:(nr+1L),2:(nc+1L)] <- replace(z,is.na(z),0)
    middle <- pad[2:(nr+1L),2:(nc+1L)]
    side <- function(neighbor,side_length) {
      hit <- middle>0 & middle!=neighbor
      if (any(hit)) perimeter <<- perimeter + tabulate(as.integer(middle[hit]),nbins=n)*side_length
    }
    side(pad[1:nr,2:(nc+1L)],rs[1])
    side(pad[3:(nr+2L),2:(nc+1L)],rs[1])
    side(pad[2:(nr+1L),1:nc],rs[2])
    side(pad[2:(nr+1L),3:(nc+2L)],rs[2])
    areas <- cnt*ca
    allp[[k]] <- data.frame(class=vals[k],patch=seq_len(n),cells=cnt,
      area=areas,perimeter=perimeter,shape_index=perimeter/(2*sqrt(pi*areas)))
    total_area <- sum(areas)
    alls[[k]] <- data.frame(class=vals[k],patches=n,area=total_area,
      mean_area=mean(areas),largest_patch_index=max(areas)/total_area,
      edge_density=sum(perimeter)/(terra::ncell(r)*ca))
  }
  list(patches=do.call(rbind,allp),summary=do.call(rbind,alls))
}

#' Compare an observed spatial statistic with a simulation ensemble
#'
#' For continuous rasters, uses nlm_variogram(); for point patterns uses
#' nlm_k(). Observed and simulated curves share requested spatial distances.
#' Envelopes are POINTWISE quantiles, not global Monte Carlo confidence bands.
#' Simulation windows and raster geometries must match observations.
#'
#' @param observed Observed raster, point result, or XY point table.
#' @param simulations List of simulated objects of the same kind.
#' @param statistic 'variogram','L','K', or 'g'.
#' @param r Spatial scales for point statistics.
#' @param lag,max_dist,angle,tolerance,pairs Variogram options.
#' @param window Observation window for point tables, or raster template.
#' @param probs Pointwise quantiles for lower and upper envelopes.
#' @param seed Random seed for sampled variogram pairs.
#' @param plot Display comparison using base graphics.
#' @return List containing a tidy curve table, observed and simulation statistics,
#'   with all scales and quantiles. Curves include 'inside' indicator.
#' @examples
#' observed <- nlm_points(method='thomas',n=100,seed=1)
#' sims <- lapply(2:21,function(s) nlm_points(method='thomas',n=100,seed=s))
#' cmp <- nlm_compare(observed,sims,statistic='L')
#' @export
nlm_compare <- function(observed,simulations,statistic=c('variogram','L','K','g'),
                        r=NULL,lag=NULL,max_dist=NULL,angle=NULL,tolerance=22.5,
                        pairs=50000,window=NULL,probs=c(0.025,0.975),
                        seed=1,plot=TRUE) {
  statistic <- match.arg(statistic)
  if (!is.list(simulations) || !length(simulations)) stop('simulations must be a nonempty list')
  if (length(probs)!=2 || any(probs<0 | probs>1) || probs[1]>=probs[2]) stop('Invalid probs')
  if (statistic=='variogram') {
    template <- .chile_diag_raster(observed)
    same <- vapply(simulations,function(s) {
      ss <- .chile_diag_raster(s)
      terra::compareGeom(template,ss,stopOnError=FALSE)
    },logical(1))
    if (!all(same)) stop('Raster geometries must match')
    fun <- function(obj,s) nlm_variogram(obj,lag=lag,max_dist=max_dist,
      angle=angle,tolerance=tolerance,pairs=pairs,seed=s)
    obs <- fun(observed,seed)
    sims <- lapply(seq_along(simulations),function(i) fun(simulations[[i]],seed+i))
    key <- paste(obs$angle,obs$distance)
    value <- 'gamma'; scale <- 'distance'
    for (z in sims) if (!identical(key,paste(z$angle,z$distance))) stop('Variogram bins differ')
  } else {
    if (is.null(window) && is.list(observed) && inherits(observed$raster,'SpatRaster')) window <- observed$raster
    if (is.null(window)) stop('Provide window for point data')
    if (inherits(window,'SpatRaster') && is.null(r)) {
      e <- terra::ext(window)
      r <- seq(min(e[2]-e[1],e[4]-e[3])/150,
               min(e[2]-e[1],e[4]-e[3])/4,length.out=30)
    }
    if (is.null(r) && is.numeric(window)) {
      r <- seq(min(diff(window[1:2]),diff(window[3:4]))/150,
               min(diff(window[1:2]),diff(window[3:4]))/4,length.out=30)
    }
    obs <- nlm_k(observed,r=r,window=window)
    sims <- lapply(simulations,function(s) nlm_k(s,r=obs$r,window=window))
    value <- statistic;scale <- 'r'
  }
  matrix_s <- vapply(sims,function(z) z[[value]],numeric(nrow(obs)))
  if (is.null(dim(matrix_s))) matrix_s <- matrix(matrix_s,ncol=1)
  lower <- apply(matrix_s,1,stats::quantile,probs=probs[1],na.rm=TRUE,names=FALSE)
  median <- apply(matrix_s,1,stats::median,na.rm=TRUE)
  upper <- apply(matrix_s,1,stats::quantile,probs=probs[2],na.rm=TRUE,names=FALSE)
  curve <- data.frame(distance=obs[[scale]],observed=obs[[value]],
    lower=lower,median=median,upper=upper,
    inside=obs[[value]]>=lower & obs[[value]]<=upper)
  if (statistic=='variogram') curve$angle <- obs$angle
  if (plot) {
    groups <- if ('angle' %in% names(curve)) unique(curve$angle) else 'all'
    oldpar <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(oldpar),add=TRUE)
    if (length(groups)>1) graphics::par(mfrow=c(ceiling(sqrt(length(groups))),ceiling(length(groups)/ceiling(sqrt(length(groups)))))
    for (g in groups) {
      d <- if (g=='all') curve else curve[curve$angle==g,,drop=FALSE]
      yr <- range(unlist(d[c('observed','lower','upper')]),finite=TRUE)
      graphics::plot(d$distance,d$observed,type='l',lwd=2,col='black',ylim=yr,
        xlab='Distance (map units)',ylab=statistic,main=if(g=='all') statistic else paste('Direction',g))
      graphics::lines(d$distance,d$lower,col='steelblue',lty=2)
      graphics::lines(d$distance,d$upper,col='steelblue',lty=2)
      graphics::lines(d$distance,d$median,col='steelblue',lwd=2)
    }
  }
  invisible(list(curve=curve,observed=obs,simulations=sims,statistic=statistic))
}
