\name{nlm_mix}
\alias{nlm_mix}
\title{Mix generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_mix(..., weights = NULL, mean = 0, sd = 1)}
\arguments{
  \item{weights}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{mean}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{sd}{See default and function-specific parameter definition; interpreted in simulation units.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
