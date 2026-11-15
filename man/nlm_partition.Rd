\name{nlm_partition}
\alias{nlm_partition}
\title{Partition generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_partition(ncol=256, nrow=256, resolution=1, n=25, method=c('voronoi','grid','bsp'), angle=0, template=NULL, seed=NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{n}{Number of requested objects or points, depending on generator.}
  \item{method}{Select the stochastic process or composition algorithm.}
  \item{angle}{Rotation in degrees counterclockwise from positive x.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
\examples{
plot(nlm_partition(ncol=80,nrow=80,n=16)$raster)
}
