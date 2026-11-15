\name{nlm_ambient}
\alias{nlm_ambient}
\title{Ambient generator}
\description{Uses spectral multiscale approximations for perlin, simplex, and value modes, not the canonical noise algorithms.}
\usage{nlm_ambient(ncol=256, nrow=256, resolution=1, extent=NULL, crs='', type='perlin', frequency=0.02, fractal=NULL, octaves=4, lacunarity=2, gain=0.5, normalize=TRUE, rescale_range=c(0,1), template=NULL, seed=NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{extent}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{crs}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{type}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{frequency}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{fractal}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{octaves}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{lacunarity}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{gain}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{normalize}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{rescale_range}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
