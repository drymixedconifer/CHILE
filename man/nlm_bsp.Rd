\name{nlm_bsp}
\alias{nlm_bsp}
\title{Bsp generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_bsp(ncol=256, nrow=256, resolution=1, n_patches=100, min_size=10, max_size=80, max_ratio=4, environment=NULL, template=NULL, overwrite=TRUE, seed=NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{n_patches}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{min_size}{Minimum geometric dimension in cells.}
  \item{max_size}{Maximum geometric dimension in cells.}
  \item{max_ratio}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{environment}{Optional raster conditioning point placement or geometry.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{overwrite}{Whether later objects replace earlier raster IDs.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
