\name{nlm_gaussianfield_matern}
\alias{nlm_gaussianfield_matern}
\title{Gaussianfield matern generator}
\description{Stationary Matérn spectral field on a periodic FFT grid; finite-domain approximation.}
\usage{nlm_gaussianfield_matern(ncol=256, nrow=256, resolution=1, range=20, sill=1, mean=0, nu=1, angle=0, range_minor=range, template=NULL, seed=NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{range}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{sill}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{mean}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{nu}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{angle}{Rotation in degrees counterclockwise from positive x.}
  \item{range_minor}{Minor correlation or cluster length in map units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
