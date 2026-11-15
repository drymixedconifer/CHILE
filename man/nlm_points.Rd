\name{nlm_points}
\alias{nlm_points}
\title{Points generator}
\description{Returns a list with point-count raster, continuous XY points, and model parameters. Fixed n conditions on the realized count. With n=NULL, applicable models generate random counts. Strauss uses approximate Metropolis relocation, inhibition uses sequential rejection.}
\usage{nlm_points(ncol=256, nrow=256, resolution=1, method=c('binomial','poisson','inhomogeneous','thomas',
                         'matern_cluster','inhibition','strauss','lgcp','grid',
                         'jittered','anisotropic_cluster','mixture'), n=100, intensity=NULL, environment=NULL, parents=12, offspring=10, cluster_sd=12, radius=15, min_dist=8, interaction=0.5, range_major=35, range_minor=8, angle=0, spacing=NULL, jitter=0.3, mixture=c('thomas','inhibition'), weights=NULL, max_attempts=100000, template=NULL, seed=NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{method}{Select the stochastic process or composition algorithm.}
  \item{n}{Number of requested objects or points, depending on generator.}
  \item{intensity}{Expected number of points per unit map area.}
  \item{environment}{Optional raster conditioning point placement or geometry.}
  \item{parents}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{offspring}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{cluster_sd}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{radius}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{min_dist}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{interaction}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{range_major}{Major correlation or cluster length in map units.}
  \item{range_minor}{Minor correlation or cluster length in map units.}
  \item{angle}{Rotation in degrees counterclockwise from positive x.}
  \item{spacing}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{jitter}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{mixture}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{weights}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{max_attempts}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{List with raster, points and parameters.}
\examples{
p <- nlm_points(ncol=80,nrow=80,n=100,method="thomas",seed=1)
plot(p$raster)
plot(p$points$x,p$points$y,asp=1)
}
