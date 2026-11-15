\name{nlm_midpoint}
\alias{nlm_midpoint}
\title{Midpoint generator}
\description{Diamond-square midpoint displacement; an enclosing 2^k+1 square is cropped to the requested output dimensions.}
\usage{nlm_midpoint(size=129, roughness=0.7, ncol=size, nrow=size, resolution=1, template=NULL, seed=NULL)}
\arguments{
  \item{size}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{roughness}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
