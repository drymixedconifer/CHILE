\name{nlm_metrics}
\alias{nlm_metrics}
\title{Connected categorical patch metrics}
\description{Summarizes patch area, perimeter, shape, largest patch dominance and edge density on a projected raster. Background NA is excluded.}
\usage{nlm_metrics(x, threshold = NULL, directions = 8)}
\arguments{\item{x}{Categorical SpatRaster or CHILE object list.}\item{threshold}{Optional continuous-field threshold to form one binary class.}\item{directions}{4 or 8 connected-cell adjacency.}}
\value{List of per-patch and per-class summary data.frames.}
\examples{b <- nlm_blob(ncol = 60, nrow = 60, n = 8, seed = 1)
nlm_metrics(b)$summary}
