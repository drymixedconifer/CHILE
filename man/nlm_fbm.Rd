\name{nlm_fbm}
\alias{nlm_fbm}
\title{Fbm generator}
\description{Spectral approximation of fractional Brownian surface, not an exact fBm simulation.}
\usage{nlm_fbm(ncol=256, nrow=256, resolution=1, H=0.7, template=NULL, seed=NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{H}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
