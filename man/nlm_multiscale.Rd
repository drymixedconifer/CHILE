\name{nlm_multiscale}
\alias{nlm_multiscale}
\title{Multiscale generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_multiscale(ncol = 256, nrow = 256, resolution = 1, ranges = c(8, 30, 90), weights = c(1, 1, 1), mean = 0, sd = 1, template = NULL, seed = NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{ranges}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{weights}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{mean}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{sd}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
