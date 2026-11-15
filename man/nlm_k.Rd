\name{nlm_k}
\alias{nlm_k}
\title{Ripley's K, centered L and annular pair correlation}
\description{Translation edge-corrected homogeneous point-pattern statistics in a rectangular projected window. Inhomogeneity can mimic clustering.}
\usage{nlm_k(x, r = NULL, window = NULL)}
\arguments{\item{x}{CHILE point result, XY data.frame or two-column XY matrix.}\item{r}{Strictly increasing distance vector in coordinate units.}\item{window}{Rectangular bounds xmin,xmax,ymin,ymax or a projected raster. Inferred from CHILE point result.}}
\value{data.frame with columns r, K, L, g and n_pairs.}
\examples{p <- nlm_points(ncol = 60, nrow = 60, method = "thomas", n = 80, seed = 1)
k <- nlm_k(p)
plot(k$r, k$L, type = "l")}
