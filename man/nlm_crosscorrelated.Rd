\name{nlm_crosscorrelated}
\alias{nlm_crosscorrelated}
\title{Crosscorrelated generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_crosscorrelated(ncol = 256, nrow = 256, resolution = 1, rho = 0.6, range1 = 30, range2 = 10, template = NULL, seed = NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{rho}{Target zero-lag Pearson correlation between output fields.}
  \item{range1}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{range2}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
