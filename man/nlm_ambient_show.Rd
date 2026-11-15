\name{nlm_ambient_show}
\alias{nlm_ambient_show}
\title{Compare procedural noise settings}
\description{Render a gallery varying up to two ambient parameters, keeping the seed constant.}
\usage{nlm_ambient_show(type="simplex",fractal="none",frequency=0.02,
  perturbation_amplitude=0,ncol=100,nrow=100,seed=42,...)}
\arguments{
\item{type}{One or more types; see \code{nlm_ambient}.}
\item{fractal}{One or more fractal families.}
\item{frequency}{One or more frequencies.}
\item{perturbation_amplitude}{One or more warp strengths; automatically selects fractal warping for positive values unless overridden in \code{...}.}
\item{ncol,nrow}{Gallery preview dimensions.}
\item{seed}{Random seed held fixed across panels.}
\item{...}{Other scalar arguments forwarded to \code{nlm_ambient}.}
}
\details{At most two of the four sweep parameters can vary; at most 24 panels. The same seed provides reproducibility but different algorithms are not necessarily paired realizations.}
\value{Invisibly returns the named list of displayed rasters.}
\examples{
\dontrun{
nlm_ambient_show(type=c("perlin","simplex","worley","cubic"),
                 fractal=c("none","fbm"))
nlm_ambient_show(frequency=c(0.01,0.03),
                 perturbation_amplitude=c(0,10,30))
}}
