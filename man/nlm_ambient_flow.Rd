\name{nlm_ambient_flow}
\alias{nlm_ambient_flow}
\title{Derive a vector field from procedural noise}
\description{Derive horizontal gradient or 2D curl vectors from an ambient noise surface using finite differences in map units.}
\usage{nlm_ambient_flow(mode="curl", ...)}
\arguments{
\item{mode}{Either \code{gradient} or \code{curl}.}
\item{...}{All arguments forwarded to \code{nlm_ambient}.}
}
\details{A gradient points toward increasing field values; curl is the 90-degree rotation of that gradient. The interior curl of a scalar field is theoretically divergence-free. This uses centered finite differences (one-sided on boundaries) and is not a hydrological flow-routing model.}
\value{Two-layer \code{SpatRaster} with layers \code{dx} and \code{dy}.}
\examples{
v <- nlm_ambient_flow(mode="curl",ncol=50,nrow=50,seed=1)
}
