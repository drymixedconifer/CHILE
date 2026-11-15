\name{nlm_radial}
\alias{nlm_radial}
\title{Radial generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_radial(ncol = 256, nrow = 256, resolution = 1, center = NULL, decay = 1, invert = FALSE, template = NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{center}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{decay}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{invert}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
