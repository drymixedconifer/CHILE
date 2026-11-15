\name{nlm_variogram}
\alias{nlm_variogram}
\title{Sampled empirical directional semivariogram}
\description{Randomly samples raster-cell pairs and computes semivariance by lag distance and angle. Use projected rasters.}
\usage{nlm_variogram(x, lag = NULL, max_dist = NULL, angle = NULL,
  tolerance = 22.5, pairs = 50000, seed = NULL)}
\arguments{\item{x}{Projected continuous SpatRaster.}\item{lag, max_dist}{Bin width and maximum distance in map units.}\item{angle}{Directions in degrees counterclockwise from x; NULL means omnidirectional.}\item{tolerance}{Directional half-width, degrees.}\item{pairs}{Sampled cell pairs.}\item{seed}{Optional seed.}}
\value{A data.frame with angle, distance, gamma and n_pairs.}
\examples{x <- nlm_anisotropic(ncol = 60, nrow = 60, seed = 1)
v <- nlm_variogram(x, angle = c(0, 90), pairs = 10000, seed = 2)
head(v)}
