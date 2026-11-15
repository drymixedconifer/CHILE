\name{nlm_show}
\alias{nlm_show}
\title{Show generator}
\description{No arguments displays a gallery of representative generators. Base graphics only.}
\usage{nlm_show(..., methods=NULL, family=NULL, ncol=90, nrow=90, seed=42)}
\arguments{
  \item{methods}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{family}{See default and function-specific parameter definition; interpreted in simulation units.}
  \item{ncol}{Number of columns when template is absent.}
  \item{nrow}{Number of rows when template is absent.}
  \item{seed}{Optional seed for reproducibility.}
}
\value{A terra SpatRaster, or a list containing a SpatRaster and associated object information, depending on the generator.}
\examples{
nlm_show(methods=c("nlm_anisotropic","nlm_mosaic","nlm_points"))
}
