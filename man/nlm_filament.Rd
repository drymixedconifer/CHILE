\name{nlm_filament}
\alias{nlm_filament}
\title{Filament generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_filament(ncol=256, nrow=256, resolution=1, n=20, length=c(30,100), width=c(1,5), mode=c('walk','straight','field','branching'), tortuosity=0.2, branching=0.02, environment=NULL, template=NULL, seed=NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{n}{Number of requested objects or points, depending on generator.}
  \item{length}{Minimum and maximum filament lengths in cells.}
  \item{width}{Minimum and maximum width in cells.}
  \item{mode}{Select the filament behavior.}
  \item{tortuosity}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{branching}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{environment}{Optional raster conditioning point placement or geometry.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
\examples{
plot(nlm_filament(ncol=80,nrow=80,n=8)$raster)
}
