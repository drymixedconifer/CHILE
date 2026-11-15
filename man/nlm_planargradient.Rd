\name{nlm_planargradient}
\alias{nlm_planargradient}
\title{Planargradient generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_planargradient(ncol=256, nrow=256, resolution=1, angle=0, template=NULL)}
\arguments{
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{resolution}{Cell resolution for new simulation rasters.}
  \item{angle}{Rotation in degrees counterclockwise from positive x.}
  \item{template}{Optional terra SpatRaster setting grid, extent, resolution and CRS.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
