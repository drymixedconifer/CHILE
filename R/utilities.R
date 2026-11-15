# Discoverable list of generators and operations.
nlm_list <- function(family=NULL) {
  catalog <- data.frame(
    function_name=c('nlm_random','nlm_gaussianfield','nlm_gaussianfield_matern',
      'nlm_anisotropic','nlm_multiscale','nlm_nonstationary','nlm_spectral',
      'nlm_fbm','nlm_midpoint','nlm_ambient','nlm_ambient_flow','nlm_planargradient',
      'nlm_radial','nlm_distance','nlm_percolation','nlm_neighborhood',
      'nlm_curdle','nlm_crosscorrelated','nlm_mosaic','nlm_points',
      'nlm_polygon','nlm_blob','nlm_filament','nlm_partition','nlm_bsp',
      'nlm_mix','nlm_transform','nlm_compose','nlm_show','nlm_ambient_show','nlm_list'),
    family=c(rep('field',19),'points',rep('geometry',5),rep('operation',3),rep('utility',3)),
    stringsAsFactors=FALSE
  )
  if (!is.null(family)) catalog <- catalog[catalog$family %in% family,,drop=FALSE]
  rownames(catalog) <- NULL;catalog
}

#' Display CHILE simulations
#'
#' Without arguments, displays a deterministic gallery of representative
#' generators. Supply raster objects, point-process results, or names of
#' generators to compare selected simulations. Base graphics only.
#'
#' @param ... Rasters, result lists, or generator names.
#' @param methods Optional vector of generator names, e.g. `c('nlm_mosaic','nlm_blob')`.
#' @param family Optional `field`, `points`, or `geometry`; used with methods.
#' @param ncol,nrow Gallery dimensions.
#' @param seed Gallery seed. The previous random state is restored.
#' @return Invisible list of displayed simulation outputs.
#' @examples
#' nlm_show(methods=c('nlm_anisotropic','nlm_mosaic','nlm_filament'))
#' @export
nlm_show <- function(...,methods=NULL,family=NULL,ncol=90,nrow=90,seed=42) {
  args <- list(...)
  if (is.null(methods) && !length(args)) {
    methods <- c('nlm_random','nlm_gaussianfield','nlm_anisotropic',
      'nlm_multiscale','nlm_nonstationary','nlm_spectral',
      'nlm_percolation','nlm_curdle','nlm_mosaic',
      'nlm_points','nlm_blob','nlm_filament')
  }
  if (!is.null(methods)) {
    if (!is.null(family)) methods <- intersect(methods,nlm_list(family)$function_name)
    oldseed <- if(exists('.Random.seed',.GlobalEnv,inherits=FALSE)) get('.Random.seed',.GlobalEnv) else NULL
    on.exit({if(is.null(oldseed)) {if(exists('.Random.seed',.GlobalEnv,inherits=FALSE)) rm('.Random.seed',envir=.GlobalEnv)}
             else assign('.Random.seed',oldseed,envir=.GlobalEnv)},add=TRUE)
    set.seed(seed)
    for (method in methods) {
      if (!method %in% nlm_list()$function_name) stop('Unknown generator: ',method)
      fun <- get(method,mode='function')
      ff <- names(formals(fun))
      params <- list()
      if ('ncol' %in% ff) params$ncol <- ncol
      if ('nrow' %in% ff) params$nrow <- nrow
      if (method=='nlm_points') {params$n<-150;params$method<-'thomas'}
      if (method=='nlm_filament') params$n<-12
      args[[length(args)+1L]] <- do.call(fun,params)
    }
    names(args)[seq.int(length(args)-length(methods)+1L,length(args))] <- methods
  }
  if (!length(args)) stop('Nothing to display')
  nm <- names(args)
  if (is.null(nm)) nm <- rep('',length(args))
  nm[nm==''] <- paste0('Landscape ',which(nm==''))
  opar <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(opar),add=TRUE)
  cols <- ceiling(sqrt(length(args)));rows <- ceiling(length(args)/cols)
  graphics::par(mfrow=c(rows,cols),mar=c(2,2,3,1))
  for(i in seq_along(args)) {
    x <- args[[i]]
    if (is.list(x) && inherits(x$raster,'SpatRaster')) {
      if (!is.null(x$points) && nrow(x$points)) {
        terra::plot(x$raster,main=nm[i],legend=FALSE)
        graphics::points(x$points$x,x$points$y,pch=16,cex=.28,col='black')
      } else terra::plot(x$raster,main=nm[i],legend=FALSE)
    } else if(inherits(x,'SpatRaster')) terra::plot(x[[1]],main=nm[i],legend=FALSE)
    else stop('Unsupported nlm_show input: ',nm[i])
  }
  invisible(args)
}
