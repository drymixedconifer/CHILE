\name{nlm_compose}
\alias{nlm_compose}
\title{Compose generator}
\description{Generate a CHILE neutral-landscape object or operate on compatible raster fields. Spatial outputs use terra. External environmental data are optional.}
\usage{nlm_compose(background, objects, method=c('replace','add','multiply'), transparent=NA_real_)}
\arguments{
  \item{background}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{objects}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{method}{Select the stochastic process or composition algorithm.}
  \item{transparent}{See default and function-specific parameter definition; interpreted in simulation units.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
\examples{
x <- nlm_random(ncol=80,nrow=80)
p <- nlm_polygon(ncol=80,nrow=80,n=5)
plot(nlm_compose(x,p))
}
