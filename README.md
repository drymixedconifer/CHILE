# CHILE v0.9 — Chaos In-Silico Landscape Engine

Standalone generators for spatial random fields, point processes, geometric objects, and composition. **terra + base R**.

## Install

```r
install.packages('terra')
install.packages('CHILE_0.9.0.tar.gz', repos=NULL, type='source')
library(CHILE)
```

## Start here

```r
nlm_list()
nlm_show()
?nlm_mosaic
?nlm_points

x <- nlm_anisotropic(range_major=60, range_minor=8, angle=45)
p <- nlm_points(method='thomas', n=200, seed=1)
nlm_show(x, p)

# Each geometry generator works by itself:
unit <- nlm_polygon(n=20, angle=c(20,60))
paths <- nlm_filament(n=15, mode='branching')

# Conditional placement and separate composition:
unit2 <- nlm_polygon(environment=x, n=20)
landscape <- nlm_compose(x,unit2)
```

## Conventions

- `ncol`, `nrow`, `resolution`, `template`, `seed` whenever stochastic and applicable.
- All angles are degrees counterclockwise from positive x. Anisotropy is specified by `range_major`, `range_minor`, and `angle` in map units.
- Fields return a `SpatRaster`; points return a list with `$raster`, `$points` (continuous XY), and `$parameters`; geometry returns `$raster` and either `$objects` or `$paths`.
- `n` in `nlm_points()` requests a fixed count. Set `n=NULL` to draw naturally random counts for processes supporting this. Clusters reaching domain edges may have fewer retained points.
- `nlm_gaussianfield_matern()` is a spectral approximation with periodic FFT boundaries. `nlm_fbm()` is a spectral approximation, **not exact fractional Brownian motion**. `nlm_midpoint()` is a multiscale approximation, not a diamond-square implementation. `nlm_ambient()` supplies multiscale spectral noise approximations for Perlin/simplex/value labels, **not their canonical noise algorithms**.
- `nlm_points(method='strauss')` uses fixed-N Metropolis updates, not exact Gibbs sampling. `method='inhibition'` uses sequential inhibition, not a stationary hard-core Gibbs process.
- Large rasters and high-density point configurations can be computationally expensive; gallery defaults are deliberately small.

## Development status

This is a package draft prepared without an available R interpreter. It has static checks and test scripts, but **has not been executed in R**. Validate with `R CMD check` before scientific use or public release.
