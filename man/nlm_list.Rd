\name{nlm_list}
\alias{nlm_list}
\title{List generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_list(family=NULL)}
\arguments{
  \item{family}{See default and function-specific parameter definition; interpreted in simulation units.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
\examples{
nlm_list()
}
