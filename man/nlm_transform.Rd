\name{nlm_transform}
\alias{nlm_transform}
\title{Transform generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_transform(nlm, distribution = c("normal", "uniform", "lognormal",
                                                  "gamma", "beta", "exponential"), params = NULL, zero_prob = 0, seed = NULL)}
\arguments{
  \item{nlm}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{distribution}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{params}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{zero_prob}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
