\name{nlm_anisotropic}
\alias{nlm_anisotropic}
\title{Anisotropic generator}
\description{Anisotropy: range_major and range_minor in grid map units; angle counterclockwise from the positive x axis.}
\usage{nlm_anisotropic(ncol = 256, nrow = 256, resolution = 1, range_major = 35, range_minor = 8, angle = 0, mean = 0, sd = 1, template = NULL, seed = NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{range_major}{Major correlation or cluster length in map units.}
  \item{range_minor}{Minor correlation or cluster length in map units.}
  \item{angle}{Rotation in degrees counterclockwise from positive x.}
  \item{mean}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{sd}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
\examples{
plot(nlm_anisotropic(ncol=80,nrow=80,range_major=30,range_minor=5,angle=45))
}
