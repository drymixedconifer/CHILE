\name{nlm_mosaic}
\alias{nlm_mosaic}
\title{Mosaic generator}
\description{Nearest-seed Voronoi mosaic, optionally in an anisotropic distance metric.}
\usage{nlm_mosaic(ncol=256, nrow=256, resolution=1, npts=20, angle=0, range_major=1, range_minor=1, template=NULL, seed=NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{npts}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{angle}{Rotation in degrees counterclockwise from positive x.}
  \item{range_major}{Major correlation or cluster length in map units.}
  \item{range_minor}{Minor correlation or cluster length in map units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
\examples{
x <- nlm_mosaic(ncol=80,nrow=80,npts=12,angle=30)
plot(x)
}
