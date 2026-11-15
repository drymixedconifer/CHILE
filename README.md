# CHILE

**Chaos in Silico Landscape Engine**  
*Experimental R package · v0.9.0*

CHILE generates synthetic landscapes from spatial random fields, point processes, and geometric objects. It is designed for controlled experiments in landscape ecology: change the structure of a landscape, hold other properties constant, and examine the consequences.

**Bring your own environment—or don't.** Every generator can run independently. An optional `terra` raster can define the simulation grid or condition the placement and development of objects. CHILE uses **`terra` and base R** for its spatial operations. No `sf`, DEM, D8 flow-direction raster, or external spatial dataset is required.

## Install and explore

From the source package archive:

```r
install.packages("CHILE_0.9.0.tar.gz", repos = NULL, type = "source")
library(CHILE)
library(terra)  # enables the familiar plot(x) workflow for SpatRaster objects

nlm_list()                          # available generators and families
nlm_show()                          # gallery of representative landscapes
nlm_show(methods = c("nlm_anisotropic", "nlm_points", "nlm_filament"))
help("nlm_points", package = "CHILE")
```

Use `terra::plot(x)` explicitly if `plot(x)` does not dispatch to terra's raster method. Individual function help pages are included as `.Rd` files; a help-viewer *internal server error* is not necessarily a failure of the generator. The installed package version and its help server should be checked separately.

## The model families

| Family | What it generates | Core return |
|---|---|---|
| **Fields** | Continuous, binary, or categorical spatial variation | `SpatRaster` |
| **Points** | Locations of discrete events or individuals | `$points`, `$raster`, `$parameters` |
| **Geometry** | Polygons, blobs, corridors, and partitions | `$raster` plus `$objects` or `$paths` |
| **Operations** | Transformations, mixtures, and overlays | `SpatRaster` |

Most generators accept `ncol`, `nrow`, `resolution`, `template`, and `seed` where relevant. Without `template`, the default grid is **256 × 256 cells**, resolution 1, and extent starting at (0, 0). With `template`, its grid geometry is used. Physical interpretations of range, distance, and angle require sensible map units—prefer a projected CRS for environmental simulations.

**Angles** are in degrees counterclockwise from positive x. Field anisotropy is specified with `range_major`, `range_minor`, and `angle`; a 1:1 range ratio is isotropic. A specified correlation scale is a model parameter, not a guarantee about any finite realization.

## 1. Spatial random fields

| Generator | Key parameters | Pattern / interpretation |
|---|---|---|
| `nlm_random()` | `minval`, `maxval` | Independent uniform noise |
| `nlm_gaussianfield()` | `range`, `range_minor`, `angle`, `sill`, `mean` | Spatially correlated Gaussian field; supports directional structure |
| `nlm_gaussianfield_matern()` | `range`, `range_minor`, `nu`, `angle`, `sill` | Matérn-like field with smoothness control |
| `nlm_anisotropic()` | `range_major`, `range_minor`, `angle`, `mean`, `sd` | Oriented Gaussian correlation |
| `nlm_multiscale()` | `ranges`, `weights`, `mean`, `sd` | Independent spatial scales combined into one field |
| `nlm_nonstationary()` | `range_small`, `range_large`, `transition`, `variance_gradient` | Spatially varying smoothness and amplitude |
| `nlm_spectral()` | `beta`, `low_cutoff`, `mean`, `sd` | Spectral texture; higher `beta` emphasizes broad features |
| `nlm_fbm()` | `H` | Approximate fractal Brownian-like surface |
| `nlm_midpoint()` | `roughness`, `size` | Multiscale terrain-like approximation |
| `nlm_ambient()` | `type`, `frequency`, `octaves`, `lacunarity`, `gain` | Procedural noise approximations and distance-based texture |
| `nlm_planargradient()` | `angle` | Linear directional gradient |
| `nlm_radial()` | `center`, `decay`, `invert` | Distance-based radial gradient, scaled to 0–1 |
| `nlm_distance()` | `n_features`, `features`, `invert` | Distance to synthetic or supplied XY feature points |
| `nlm_percolation()` | `prob` | Independent binary occupancy |
| `nlm_neighborhood()` | `beta`, `intercept`, `sweeps`, `neighbors` | Binary neighborhood-dependent field (Ising-like Gibbs updates) |
| `nlm_curdle()` | `levels`, `p` | Recursive hierarchical binary pattern |
| `nlm_crosscorrelated()` | `rho`, `range1`, `range2` | Two correlated Gaussian layers |
| `nlm_mosaic()` | `npts`, `range_major`, `range_minor`, `angle` | Nearest-seed (Voronoi-like) categorical mosaic |

`nlm_mosaiccluster()` remains an alias for `nlm_mosaic()`; its historical `cluster_sd` argument does not alter the result.

**Start with anisotropy and scale:**

```r
x <- nlm_anisotropic(range_major = 60, range_minor = 8,
                     angle = 45, seed = 1)
y <- nlm_multiscale(ranges = c(5, 25, 90),
                    weights = c(1, 0.7, 1.5), seed = 2)
nlm_show(x, y)
```

**Joint environmental variation:**

```r
climate <- nlm_crosscorrelated(rho = -0.7, range1 = 40,
                                range2 = 15, seed = 1)
terra::plot(climate)
```

## 2. Spatial point processes

`nlm_points()` generates **continuous XY locations**, retaining subcell precision; `$raster` contains counts per cell. Set `method` to choose the process.

| `method` | Main controls | Behavior |
|---|---|---|
| `"binomial"` | `n` | Fixed-count complete spatial randomness |
| `"poisson"` | `n = NULL`, `intensity` | Random count with independently located points |
| `"inhomogeneous"` | `environment`, `n` or `intensity` | Location intensity weighted by a nonnegative spatial surface |
| `"thomas"` | `parents`, `offspring`, `cluster_sd` | Gaussian offspring clusters |
| `"matern_cluster"` | `parents`, `offspring`, `radius` | Circular offspring clusters |
| `"anisotropic_cluster"` | `parents`, `range_major`, `range_minor`, `angle` | Oriented Gaussian clusters |
| `"inhibition"` | `min_dist`, `max_attempts` | Sequential minimum-distance inhibition |
| `"strauss"` | `radius`, `interaction`, `max_attempts` | Approximate fixed-count Gibbs interaction; 0–1 repulsion |
| `"lgcp"` | `range_major`, `range_minor`, `angle`, `intensity` | Cox-like process driven by a latent spatial field |
| `"grid"` | `spacing`, `angle` | Regularly spaced locations |
| `"jittered"` | `spacing`, `angle`, `jitter` | Disturbed regular grid |
| `"mixture"` | `mixture`, `weights` | Mixture of component process families |

Other shared arguments include `n`, `environment`, `seed`, and `template`. `n` defaults to 100; it is **a requested fixed count**, not an intensity. To simulate an unconditional Poisson count, use `n = NULL` and `intensity` (per squared map unit). Inhibition, boundary clipping, and some regular-grid configurations can yield a different realized count; inspect `p$parameters$n_realized`.

```r
# Compare aggregation, repulsion, and anisotropic clustering
p1 <- nlm_points(method = "thomas", n = 250,
                 parents = 12, cluster_sd = 10, seed = 1)
p2 <- nlm_points(method = "inhibition", n = 250,
                 min_dist = 6, seed = 1)
p3 <- nlm_points(method = "anisotropic_cluster", n = 250,
                 range_major = 30, range_minor = 5,
                 angle = 35, seed = 1)
nlm_show(p1, p2, p3)

head(p1$points)       # id, x, y
p1$parameters         # requested and realized count
terra::plot(p1$raster) # cell counts
plot(p1$points$x, p1$points$y, asp = 1, pch = 16)
```

**Condition point placement on a synthetic or observed surface:**

```r
env <- nlm_anisotropic(range_major = 45, range_minor = 10, seed = 3)
w <- nlm_transform(env, distribution = "uniform")
p <- nlm_points(method = "inhomogeneous", n = 400,
                environment = w, seed = 4)
nlm_show(w, p)
```

The `environment` supplied as point-location weights must be **nonnegative**. In the current implementation, `inhomogeneous` with `n = NULL` is a spatially weighted random-count placement model, while `lgcp` uses raster-cell random intensities. These are useful experimental constructions, but this development build should be validated before making strict distributional claims. The Strauss and sequential-inhibition methods are approximations, not exact stationary Gibbs samplers.

## 3. Geometric landscape objects

| Generator | Main controls | Outputs |
|---|---|---|
| `nlm_polygon()` | `n`, `min_size`, `max_size`, `max_ratio`, `angle`, `skew`, `taper` | Rotated straight-edged polygons; `$objects` |
| `nlm_blob()` | `n`, `min_size`, `max_size`, `max_ratio`, `irregularity`, `complexity` | Irregular elliptical patches; `$objects` |
| `nlm_filament()` | `n`, `length`, `width`, `mode`, `tortuosity`, `branching` | Linear corridors and branching networks; `$paths` |
| `nlm_partition()` | `n`, `method`, `angle` | Voronoi, rotated grid, or recursive rectangular partition |
| `nlm_bsp()` | `n_patches`, `min_size`, `max_size`, `max_ratio` | Convenience wrapper for axis-aligned rectangle placement |

Filament modes are `"straight"`, `"walk"`, `"branching"`, and `"field"`. `"field"` generates a synthetic surface internally if none is supplied; it follows locally favorable cells, **not** a hydrologically routed stream network. Polygon and filament widths and lengths are specified in cells, unlike field correlation ranges on a template, which follow map units. Existing geometry outputs retain object IDs; overlap may make some objects partially or entirely invisible in the final raster.

```r
polys <- nlm_polygon(n = 30, angle = c(15, 65), seed = 1)
blobs <- nlm_blob(n = 30, irregularity = 0.45, seed = 2)
roads <- nlm_filament(n = 12, mode = "walk",
                      tortuosity = 0.12, seed = 3)
units <- nlm_partition(n = 30, method = "voronoi", seed = 4)
nlm_show(polys, blobs, roads, units)
```

For environmentally biased **placement**, supply `environment` to `nlm_polygon()` or `nlm_blob()`; this is distinct from conditioning the shape of individual patches. In `nlm_filament()`, `environment` also affects starting locations and, in `"field"` mode, the movement field.

## 4. Transform, mix, compose

| Function | Main arguments | Purpose |
|---|---|---|
| `nlm_transform()` | `distribution`, `params`, `zero_prob` | Rank-based marginal transformation; normal, uniform, lognormal, gamma, beta, exponential |
| `nlm_mix()` | `...`, `weights`, `mean`, `sd` | Standardize and combine continuous fields |
| `nlm_compose()` | `background`, `objects`, `method`, `transparent` | Overlay on a shared grid using `"replace"`, `"add"`, or `"multiply"` |

```r
set.seed(42)
a <- nlm_anisotropic(range_major = 50, range_minor = 10, angle = 25)
b <- nlm_spectral(beta = 3)
env <- nlm_mix(a, b, weights = c(0.7, 0.3))
objects <- nlm_polygon(template = env, n = 25, min_size = 12,
                       max_size = 35, seed = 7)
landscape <- nlm_compose(env, objects, method = "replace")
nlm_show(env, objects, landscape)
```

`nlm_transform()` preserves ranks except where zero inflation adds ties; it does not generally preserve Pearson spatial covariance. `nlm_compose()` operates on **values**: when replacing, polygon IDs replace underlying field values. It does not automatically apply ecological meaning or interpret object IDs as treatment intensities.

## 5. See and document what you're generating

```r
nlm_list()                               # all available entries
nlm_list(family = "points")              # point-process family
nlm_show()                               # full default gallery
nlm_show(methods = c("nlm_mosaic", "nlm_blob"))

# Compare objects generated with different parameters
low  <- nlm_anisotropic(range_major = 20, range_minor = 20, seed = 10)
high <- nlm_anisotropic(range_major = 70, range_minor = 5,
                        angle = 45, seed = 10)
nlm_show(isotropic = low, anisotropic = high)

?nlm_mosaic
?nlm_points
?nlm_show
```

`nlm_show()` uses base graphics and `terra::plot()`, with the default gallery drawn at 90 × 90 cells for speed. It returns the displayed objects invisibly. Use `nlm_show(methods = ...)` to generate named methods with defaults; to compare parameter settings, first generate the objects and pass them to `nlm_show()`.

## Scientific and development status

CHILE is a **neutral/synthetic landscape experiment engine**, not a model of ecological mechanisms by default. Statistical controls specify generating processes; realized patch counts, connectivity, spatial covariance, and anisotropy should be measured when those properties are the response of interest.

Several methods intentionally approximate classical processes: Matérn generation uses a spectral approximation with periodic boundaries; `nlm_fbm()` is not exact fractional Brownian motion; `nlm_midpoint()` is not a diamond-square implementation; and `nlm_ambient(type = "perlin" | "simplex" | "value")` uses multiscale approximations rather than canonical noise algorithms. These are candidates for validation or replacement before a stable release.

**v0.9.0 is a development draft.** It has not been run through a verified R runtime or `R CMD check` in the build environment. Treat results as exploratory until package checks, method-level tests, and spatial diagnostics are completed.
