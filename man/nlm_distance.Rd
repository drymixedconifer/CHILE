\name{nlm_distance}
\alias{nlm_distance}
\title{Distance generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_distance(ncol = 256, nrow = 256, resolution = 1, n_features = 10, features = NULL, invert = FALSE, template = NULL, seed = NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{n_features}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{features}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{invert}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
