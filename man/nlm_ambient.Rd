\name{nlm_ambient}
\alias{nlm_ambient}
\title{Native procedural noise with ambient}
\description{Generate independent synthetic fields using the native \pkg{ambient} noise and pattern algorithms, with fractal combinations, domain warping, and coordinate-based anisotropy. No environmental input is required.}
\usage{
nlm_ambient(ncol=256, nrow=256, resolution=1, extent=NULL, crs="",
  type="simplex", frequency=0.02, fractal=NULL, octaves=4,
  lacunarity=2, gain=0.5, perturbation="none",
  perturbation_amplitude=1, angle=0, anisotropy=1,
  interpolator="quintic", distance="euclidean", value="cell",
  distance_ind=c(1,2), jitter=0.45,
  blue_sd=10, blue_seed_frac=0.1,
  normalize=TRUE, rescale_range=c(0,1), template=NULL, seed=NULL)
}
\arguments{
\item{ncol,nrow}{Grid columns and rows (default 256 each).}
\item{resolution}{Cell size; defaults to one.}
\item{extent,crs}{Optional raster extent and CRS for new grids.}
\item{type}{One of \code{perlin}, \code{simplex}, \code{cubic}, \code{value}, \code{white}, \code{worley}, \code{blue}, \code{waves}, \code{spheres}, \code{checkerboard}.}
\item{frequency}{Noise frequency in inverse \emph{raster cells}; 0.02 corresponds to roughly 50 cells per basic wavelength.}
\item{fractal}{\code{NULL} or \code{none} (one octave), \code{fbm}, \code{billow}, or \code{rigid-multi}. The previous \code{rigidmulti} spelling is accepted.}
\item{octaves}{Number of combined layers when fractal mode is enabled.}
\item{lacunarity}{Frequency multiplier between octaves.}
\item{gain}{Relative amplitude multiplier between octaves.}
\item{perturbation}{\code{none}, \code{normal}, or \code{fractal}; use \code{fractal} for multiscale domain warping.}
\item{perturbation_amplitude}{Magnitude of coordinate displacement, in cells.}
\item{angle}{Major anisotropy axis, degrees counterclockwise from the positive x direction.}
\item{anisotropy}{Major-to-minor spatial-scale ratio, at least one; 1 means isotropic.}
\item{interpolator}{\code{linear}, \code{hermite}, \code{quintic} for Perlin and value noise.}
\item{distance}{For Worley: \code{euclidean}, \code{manhattan}, or \code{natural}.}
\item{value}{For Worley: \code{cell}, \code{distance}, \code{distance2}, \code{distance2add}, \code{distance2sub}, \code{distance2mul}, \code{distance2div}.}
\item{distance_ind}{Two neighbor ranks for Worley distance2 variants.}
\item{jitter}{Worley point displacement, default 0.45.}
\item{blue_sd,blue_seed_frac}{Blue-noise generation parameters; blue is limited to 65536 cells because it is expensive.}
\item{normalize}{Rescale observed minimum and maximum to \code{rescale_range}.}
\item{rescale_range}{Two values giving the rescaled output range.}
\item{template}{Optional existing \code{terra::SpatRaster} defining dimensions and location.}
\item{seed}{Optional integer seed.}
}
\details{The native \pkg{ambient} \code{gen_*} functions evaluate noise on CHILE's cell coordinates. CHILE applies coordinate warping before an affine rotation/stretch and builds fractal combinations explicitly. Fractal formulas follow conventional procedural definitions but should not be interpreted as exact fractional Brownian random fields. The normalized realization need not have a specified covariance. \code{blue} ignores the frequency parameter and does not support domain warping, fractal synthesis, or anisotropy; white and geometric patterns also exclude fractal synthesis.}
\value{A single-layer \code{terra::SpatRaster}.}
\examples{
x <- nlm_ambient(type="simplex", frequency=0.02, seed=1)
y <- nlm_ambient(type="worley", value="distance2sub", seed=1)
z <- nlm_ambient(type="perlin", fractal="billow", angle=45,
                 anisotropy=5, perturbation="fractal",
                 perturbation_amplitude=20, seed=1)
}
