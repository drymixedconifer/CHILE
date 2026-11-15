\name{nlm_neighborhood}
\alias{nlm_neighborhood}
\title{Neighborhood generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_neighborhood(ncol = 256, nrow = 256, resolution = 1, beta = 0.7, intercept = 0, sweeps = 25, neighbors = 4, initial_prob = 0.5, template = NULL, seed = NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{beta}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{intercept}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{sweeps}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{neighbors}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{initial_prob}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
