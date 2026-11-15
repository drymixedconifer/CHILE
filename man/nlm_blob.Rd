\name{nlm_blob}
\alias{nlm_blob}
\title{Blob generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_blob(ncol=256, nrow=256, resolution=1, n=30, min_size=8, max_size=50, max_ratio=4, irregularity=0.3, complexity=5, environment=NULL, template=NULL, overwrite=TRUE, seed=NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{n}{Number of requested objects or points, depending on generator.}
  \item{min_size}{Minimum geometric dimension in cells.}
  \item{max_size}{Maximum geometric dimension in cells.}
  \item{max_ratio}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{irregularity}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{complexity}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{environment}{Optional raster conditioning point placement or geometry.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{overwrite}{Whether later objects replace earlier raster IDs.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
\examples{
plot(nlm_blob(ncol=80,nrow=80,n=10)$raster)
}
