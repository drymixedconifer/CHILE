\name{nlm_nonstationary}
\alias{nlm_nonstationary}
\title{Nonstationary generator}
\description{Local correlation scales are blended spatially. Not an exact locally stationary covariance model.}
\usage{nlm_nonstationary(ncol = 256, nrow = 256, resolution = 1, range_small = 5, range_large = 55, transition = NULL, variance_gradient = 0, mean = 0, sd = 1, template = NULL, seed = NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{range_small}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{range_large}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{transition}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{variance_gradient}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{mean}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{sd}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
