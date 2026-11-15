# Minimal raster-native geometry. Each generator works without a supplied NLM.
.chile_centercells <- function(r,n,environment=NULL) {
  if (is.null(environment)) return(sample.int(terra::ncell(r),n,replace=TRUE))
  if (!terra::compareGeom(r,environment,stopOnError=FALSE)) stop('environment grid mismatch')
  w <- terra::values(environment)[,1]; w[!is.finite(w)] <- 0
  w <- w-min(w)
  if (!any(w>0)) w[] <- 1
  sample.int(length(w),n,replace=TRUE,prob=w)
}
.chile_geom_raster <- function(r,m,label='object_id') {
  terra::values(r) <- as.vector(t(m)); names(r) <- label; r
}

nlm_polygon <- function(ncol=256,nrow=256,resolution=1,n=30,min_size=8,max_size=50,
                        max_ratio=4,angle=c(0,180),skew=0,taper=0,
                        environment=NULL,template=NULL,overwrite=TRUE,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(template) && !is.null(environment)) template <- environment
  r <- .chile_raster(ncol,nrow,resolution,template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  m <- matrix(NA_integer_,nr,nc)
  ctr <- terra::rowColFromCell(r,.chile_centercells(r,n,environment))
  pars <- data.frame(id=seq_len(n),row=ctr[,1],col=ctr[,2],
    width=sample(seq.int(min_size,max_size),n,TRUE),ratio=runif(n,1,max_ratio),
    angle=runif(n,min(angle),max(angle)))
  for (i in seq_len(n)) {
    p <- pars[i,]; w <- p$width; h <- w/p$ratio
    x <- c(-1,1,1+skew,-1+skew)*w/2
    y <- c(-1,-1,1,1)*h/2
    x[3:4] <- x[3:4]*(1-taper)
    a <- -p$angle*pi/180
    xx <- x*cos(a)-y*sin(a)+p$col
    yy <- x*sin(a)+y*cos(a)+p$row
    cc <- seq.int(max(1,floor(min(xx))),min(nc,ceiling(max(xx))))
    rr <- seq.int(max(1,floor(min(yy))),min(nr,ceiling(max(yy))))
    if (!length(cc)||!length(rr)) next
    g <- expand.grid(row=rr,col=cc); inside <- rep(FALSE,nrow(g))
    j <- length(xx)
    for(k in seq_along(xx)) {
      hit <- ((yy[k]>g$row)!=(yy[j]>g$row)) &
        (g$col<(xx[j]-xx[k])*(g$row-yy[k])/(yy[j]-yy[k]+1e-12)+xx[k])
      inside <- xor(inside,hit); j <- k
    }
    cells <- as.matrix(g[inside,c('row','col'),drop=FALSE])
    if (!nrow(cells)) next
    if (overwrite) m[cells] <- i else {
      empty <- is.na(m[cells]); m[cells[empty,,drop=FALSE]] <- i
    }
  }
  list(raster=.chile_geom_raster(r,m),objects=pars)
}

nlm_blob <- function(ncol=256,nrow=256,resolution=1,n=30,min_size=8,max_size=50,
                     max_ratio=4,irregularity=0.3,complexity=5,
                     environment=NULL,template=NULL,overwrite=TRUE,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(template) && !is.null(environment)) template <- environment
  r <- .chile_raster(ncol,nrow,resolution,template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  m <- matrix(NA_integer_,nr,nc)
  ctr <- terra::rowColFromCell(r,.chile_centercells(r,n,environment))
  pars <- data.frame(id=seq_len(n),row=ctr[,1],col=ctr[,2],
    size=sample(seq.int(min_size,max_size),n,TRUE),ratio=runif(n,1,max_ratio),
    angle=runif(n,0,pi))
  for (i in seq_len(n)) {
    p <- pars[i,]; a <- p$size/2; b <- a/p$ratio
    radius <- ceiling(a*(1+irregularity*sum(1/seq_len(complexity))))
    rr <- seq.int(max(1,p$row-radius),min(nr,p$row+radius))
    cc <- seq.int(max(1,p$col-radius),min(nc,p$col+radius))
    if (!length(rr)||!length(cc)) next
    g <- expand.grid(row=rr,col=cc)
    x <- g$col-p$col; y <- g$row-p$row
    a <- -p$angle
    xr <- x*cos(a)+y*sin(a)
    yr <- -x*sin(a)+y*cos(a)
    theta <- atan2(yr/b,xr/a); d <- sqrt((xr/a)^2+(yr/b)^2)
    noise <- numeric(length(theta))
    for (k in seq_len(complexity)) noise <- noise+sin(k*theta+runif(1,0,2*pi))/k
    inside <- d <= pmax(0.1,1+irregularity*noise)
    cells <- as.matrix(g[inside,c('row','col'),drop=FALSE])
    if (!nrow(cells)) next
    if (overwrite) m[cells] <- i else {
      empty <- is.na(m[cells]);m[cells[empty,,drop=FALSE]] <- i
    }
  }
  list(raster=.chile_geom_raster(r,m),objects=pars)
}

nlm_filament <- function(ncol=256,nrow=256,resolution=1,n=20,length=c(30,100),
                         width=c(1,5),mode=c('walk','straight','field','branching'),
                         tortuosity=0.2,branching=0.02,environment=NULL,
                         template=NULL,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  mode <- match.arg(mode)
  if (is.null(template) && !is.null(environment)) template <- environment
  r <- .chile_raster(ncol,nrow,resolution,template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  m <- matrix(NA_integer_,nr,nc)
  centers <- terra::rowColFromCell(r,.chile_centercells(r,n,environment))
  if (mode=='field' && is.null(environment)) environment <- nlm_multiscale(template=r)
  if (!is.null(environment)) env <- matrix(terra::values(environment)[,1],nr,nc,byrow=TRUE)
  paths <- vector('list',n)
  for (i in seq_len(n)) {
    steps <- sample(seq.int(length[1],length[2]),1)
    w <- sample(seq.int(width[1],width[2]),1)
    queue <- list(list(pos=as.numeric(centers[i,]),heading=runif(1,0,2*pi),steps=steps))
    segments <- list(); b <- 1L
    while (b<=length(queue) && b<=100L) {
      pos <- queue[[b]]$pos; heading <- queue[[b]]$heading
      path <- matrix(numeric(0),ncol=2)
      for (s in seq_len(queue[[b]]$steps)) {
        rr <- round(pos[1]); cc <- round(pos[2])
        if (rr<1||rr>nr||cc<1||cc>nc) break
        path <- rbind(path,c(rr,cc))
        if (mode=='branching' && runif(1)<branching && length(queue)<100L) {
          queue[[length(queue)+1L]] <- list(pos=pos,
            heading=heading+sample(c(-1,1),1)*runif(1,0.5,1.5),
            steps=max(1L,round((queue[[b]]$steps-s)*0.6)))
        }
        if (mode %in% c('walk','branching')) heading <- heading+rnorm(1,0,tortuosity)
        if (mode=='field') {
          neighbors <- expand.grid(dr=-1:1,dc=-1:1)
          candidates <- cbind(rr+neighbors$dr,cc+neighbors$dc)
          ok <- candidates[,1]>=1 & candidates[,1]<=nr & candidates[,2]>=1 & candidates[,2]<=nc
          candidates <- candidates[ok,,drop=FALSE]
          vals <- env[candidates]
          if (!any(is.finite(vals))) break
          j <- which.min(replace(vals,!is.finite(vals),Inf))
          if (all(candidates[j,]==c(rr,cc))) break
          pos <- as.numeric(candidates[j,])
        } else pos <- pos+c(sin(heading),cos(heading))
      }
      segments[[length(segments)+1L]] <- path
      b <- b+1L
    }
    paths[[i]] <- segments
    rad <- w/2
    offsets <- expand.grid(dr=-ceiling(rad):ceiling(rad),dc=-ceiling(rad):ceiling(rad))
    offsets <- offsets[offsets$dr^2+offsets$dc^2<=rad^2,,drop=FALSE]
    for(path in segments) for(j in seq_len(nrow(path))) {
      rr <- path[j,1]+offsets$dr;cc <- path[j,2]+offsets$dc
      ok <- rr>=1 & rr<=nr & cc>=1 & cc<=nc
      cells <- cbind(rr[ok],cc[ok]);m[cells] <- i
    }
  }
  list(raster=.chile_geom_raster(r,m,'filament_id'),paths=paths)
}

nlm_partition <- function(ncol=256,nrow=256,resolution=1,n=25,
                          method=c('voronoi','grid','bsp'),angle=0,
                          template=NULL,seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  method <- match.arg(method)
  r <- .chile_raster(ncol,nrow,resolution,template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  if (method=='voronoi') return(list(raster=nlm_mosaic(npts=n,template=r)))
  if (method == "grid") {

  rr <- row(matrix(0, nr, nc))
  cc <- col(matrix(0, nr, nc))

  a <- angle * pi / 180

  x <- cc - (nc + 1) / 2
  y <- (nr + 1) / 2 - rr

  u <- x * cos(a) + y * sin(a)
  v <- -x * sin(a) + y * cos(a)

  s <- sqrt(nr * nc / n)

  gx <- floor((u - min(u)) / s)
  gy <- floor((v - min(v)) / s)

  m <- gx + (max(gx) + 1) * gy + 1

} else {
    m <- matrix(NA_integer_,nr,nc)
    parts <- list(c(1L,nr,1L,nc))
    while(length(parts)<n) {
      areas <- vapply(parts,function(p)(p[2]-p[1]+1)*(p[4]-p[3]+1),numeric(1))
      j <- which.max(areas);p <- parts[[j]]
      h <- p[2]-p[1]+1;w <- p[4]-p[3]+1
      if(h<2 && w<2) break
      byrow <- if(h<2) FALSE else if(w<2) TRUE else h>=w
      if(byrow) {
        cut <- sample(seq.int(p[1],p[2]-1),1)
        a <- c(p[1],cut,p[3],p[4]);b <- c(cut+1,p[2],p[3],p[4])
      } else {
        cut <- sample(seq.int(p[3],p[4]-1),1)
        a <- c(p[1],p[2],p[3],cut);b <- c(p[1],p[2],cut+1,p[4])
      }
      parts[[j]] <- a;parts[[length(parts)+1L]] <- b
    }
    for(i in seq_along(parts)) {p<-parts[[i]];m[p[1]:p[2],p[3]:p[4]]<-i}
  }
  list(raster=.chile_geom_raster(r,m,'partition_id'))
}

# Overlay composition; zeros/NA in objects are transparent by default.
nlm_compose <- function(background,objects,method=c('replace','add','multiply'),
                        transparent=NA_real_) {
  method <- match.arg(method)
  obj <- if (is.list(objects) && inherits(objects$raster,'SpatRaster')) objects$raster else objects
  if (!terra::compareGeom(background,obj,stopOnError=FALSE)) stop('Grid mismatch')
  r <- background[[1]]
  a <- terra::values(r)[,1];b <- terra::values(obj)[,1]
  use <- is.finite(b)
  if (!is.na(transparent)) use <- use & b!=transparent
  if (method=='replace') a[use] <- b[use]
  if (method=='add') a[use] <- a[use]+b[use]
  if (method=='multiply') a[use] <- a[use]*b[use]
  terra::values(r) <- a;names(r)<-'composed';r
}

# Original nlm_bsp concept, revised: generated independent of NLM.
nlm_bsp <- function(ncol=256,nrow=256,resolution=1,n_patches=100,
                    min_size=10,max_size=80,max_ratio=4,
                    environment=NULL,template=NULL,overwrite=TRUE,seed=NULL) {
  nlm_polygon(ncol,nrow,resolution,n=n_patches,min_size=min_size,max_size=max_size,
    max_ratio=max_ratio,angle=c(0,0),environment=environment,
    template=template,overwrite=overwrite,seed=seed)
}
