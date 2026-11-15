# CHILE

**Chaos in Silico Landscape Engine** · v0.9.2 (development)

CHILE generates synthetic landscapes using **random fields**, **point processes**, and **geometric objects**. Simulations run with no external spatial inputs; an optional `terra::SpatRaster` can establish the grid or condition spatial placement. Generate components independently, then combine them when an experiment calls for it.

> **Documentation:** [Complete parameter reference](docs/reference.md) · [Point-process methods](docs/reference.md#2-point-processes-nlm_points) · [Procedural noise / ambient](docs/reference.md#nlm_ambient--native-procedural-noise-via-ambient) · [Limitations](docs/reference.md#7-current-development-limitations)

## Install and start

```r
install.packages("CHILE_0.9.2.tar.gz", repos = NULL, type = "source")
library(CHILE)

nlm_list()    # generator catalog
nlm_show()    # gallery of representative examples

x <- nlm_anisotropic(range_major = 60, range_minor = 8,
                     angle = 35, seed = 1)
terra::plot(x)
```

CHILE requires R ≥ 4.1, `terra`, and `ambient` (plus standard R libraries). **No `sf`, DEM, D8, or cost raster is required.** A source installation needs a working R build toolchain where applicable.

## Choose a generator

| Goal | Functions | Main controls |
|---|---|---|
| Spatially correlated continuous field | `nlm_gaussianfield`, `nlm_gaussianfield_matern`, `nlm_anisotropic` | `range`, `range_major`, `range_minor`, `angle`, `nu` |
| Multiscale or nonstationary structure | `nlm_multiscale`, `nlm_nonstationary`, `nlm_spectral`, `nlm_fbm`, `nlm_midpoint` | `ranges`, `weights`, `beta`, `H`, `roughness` |
| Procedural noise and warping | `nlm_ambient`, `nlm_ambient_flow`, `nlm_ambient_show` | `type`, `frequency`, `fractal`, `octaves`, `perturbation`, `anisotropy` |
| Basic or structured gradients | `nlm_random`, `nlm_planargradient`, `nlm_radial`, `nlm_distance` | amplitude, direction, centers and distance |
| Categorical and mosaic fields | `nlm_percolation`, `nlm_neighborhood`, `nlm_curdle`, `nlm_mosaic` | occupancy, neighborhood interaction, partition geometry |
| Individuals, events, stems | `nlm_points` | `method`, `n`, `intensity`, `parents`, `cluster_sd`, `min_dist` |
| Patches and networks | `nlm_polygon`, `nlm_blob`, `nlm_filament`, `nlm_partition`, `nlm_bsp` | size, elongation, irregularity, tortuosity, branching |
| Combine or transform | `nlm_transform`, `nlm_mix`, `nlm_compose` | distribution, weights, overlay operation |

For **all signatures, default values, method choices, parameter meanings, and examples**, use the [full function reference](docs/reference.md). Some functions are compatibility aliases; see their reference entries.

## The workflow: fields → objects → composition

These operations are independent. A field is optional when generating points or geometry.

```r
# 1. Synthetic environmental heterogeneity
env <- nlm_ambient(ncol = 250, nrow = 250,
                   type = "simplex", fractal = "fbm",
                   octaves = 5, frequency = 0.02,
                   anisotropy = 3, angle = 30, seed = 10)

# 2. Spatially clustered points on the same grid
stems <- nlm_points(method = "thomas", n = 400,
                    parents = 20, cluster_sd = 12,
                    template = env, seed = 11)

# 3. Geometric patches, independent of environmental conditioning
patches <- nlm_blob(n = 35, template = env, seed = 12)

# 4. Combine explicitly
combined <- nlm_compose(env, patches, method = "replace")

nlm_show(environment = env, stem_counts = stems$raster,
         patch_ids = patches$raster, combined = combined)
```

`template` supplies grid geometry; `environment` supplies spatial weighting when a generator supports conditioning. They are **not interchangeable**. `nlm_compose()` overlays raster cell values; polygon IDs are not automatically binary masks.

## Explore procedural noise

The native `ambient` backend includes Perlin, simplex, Worley, cubic, value, white and blue noise, plus procedural patterns. Fractal synthesis and coordinate warping greatly expand the range of textures.

```r
nlm_ambient_show(
  type = c("perlin", "simplex", "worley", "cubic"),
  fractal = c("none", "fbm", "billow")
)

x <- nlm_ambient(
  type = "simplex", frequency = 0.015,
  fractal = "rigid-multi", octaves = 6,
  perturbation = "fractal", perturbation_amplitude = 25,
  anisotropy = 5, angle = 35, seed = 42
)
terra::plot(x)
```

Not every fractal/warping option applies to every native noise family; the [ambient reference](docs/reference.md#nlm_ambient--native-procedural-noise-via-ambient) specifies the combinations and exceptions.

## Explore point processes

`nlm_points()` provides fixed-count CSR (`binomial`), Poisson, inhomogeneous Poisson, Thomas and Matérn clusters, anisotropic clusters, inhibition, approximate Strauss, log-Gaussian Cox, regular/jittered grids, and mixtures. It returns **continuous XY locations** and a count raster.

```r
clustered <- nlm_points(method = "thomas", n = 500,
                        parents = 15, cluster_sd = 10, seed = 2)
spaced <- nlm_points(method = "inhibition", n = 200,
                     min_dist = 7, seed = 2)

nlm_show(clustered = clustered$raster,
         spaced = spaced$raster)
head(clustered$points)
```

`n` may specify a target, a fixed count, or a count-conditioned process depending on method. For a true Poisson realization, use `n = NULL` and an appropriate `intensity`. See [the detailed point-process reference](docs/reference.md#2-point-processes-nlm_points) before comparing modeled intensities or clustering.

## Conventions

- **Output:** fields return `SpatRaster`; points return `list(raster, points, parameters)`; object generators return lists containing `$raster` and geometry metadata or paths.
- **Spatial units:** field ranges and point-process distances use the grid's map units; polygon dimensions and filament widths/lengths are measured in cells. Use a projected CRS if units should mean meters.
- **Anisotropy:** `angle` is conventionally counterclockwise from positive x for fields and point clusters. Check geometric orientation separately because some geometry uses raster row/column coordinates.
- **Reproducibility:** supply `seed` when comparing processes. Parameter settings describe generating rules, not necessarily the exact measured properties of a finite realization.
- **Experimental status:** several named processes are finite-domain or approximate implementations. See [limitations and proposed validation](docs/reference.md#7-current-development-limitations).

## Visualization and discovery

```r
nlm_list()                        # all functions
nlm_list(family = "points")       # relevant family
nlm_show()                        # default gallery
nlm_show(methods = c("nlm_anisotropic", "nlm_points", "nlm_blob"))

horizontal <- nlm_anisotropic(angle = 0, seed = 1)
vertical   <- nlm_anisotropic(angle = 90, seed = 1)
nlm_show(horizontal = horizontal, vertical = vertical)
```

`nlm_show()` accepts generated rasters or lists with a `$raster` element. The [reference](docs/reference.md#5-discover-and-visualize-nlm_list-and-nlm_show) describes its full arguments.

## R help-pane error: diagnose documentation vs viewer

The package archive contains `.Rd` help topics, but **RStudio's “internal server error” does not establish whether the problem is the installed documentation or the Help pane's local HTTP viewer**. Try the following in a fresh R session:

```r
library(CHILE)
h <- utils::help("nlm_points", package = "CHILE")
print(h)                            # can invoke IDE help viewer
length(h)                           # should be 1 if topic is indexed
rd <- utils:::.getHelpFile(h)
cat(paste(capture.output(tools::Rd2txt(rd)), collapse = "\n"))
```

If `Rd2txt()` succeeds but the Help pane fails, the `.Rd` topic is readable and the problem is likely in the IDE's help rendering/server or installation cache. Restart R/RStudio and reinstall the package. If `length(h) == 0` or `.getHelpFile()` fails, investigate the installed help index/package archive instead. The R console error text and RStudio version can narrow it further. Avoid interpreting an IDE server failure as a broken generator.

## Reference and validation

The [complete parameter reference](docs/reference.md) retains every generator signature, default, method, usage example, and scientific qualification from the previous detailed README. The engine is a **development build**: synthetic-looking structure is not evidence that empirical covariance, clustering, anisotropy, or power spectra match target parameters. Those properties need diagnostic validation before scientific inference.

## Analyze and compare simulated landscapes (v0.9.3)

Four small diagnostics provide **scale-dependent summaries** and observed-versus-simulated comparisons. No new dependencies; `terra` and base R only. For spatial distance analysis, use projected coordinates with linear map units, not longitude/latitude.

| Function | Input | Returns | Key parameters |
|---|---|---|---|
| `nlm_variogram()` | Continuous raster | Directional empirical semivariogram (`distance`, `gamma`, `n_pairs`) | `lag`, `max_dist`, `angle`, `tolerance`, `pairs` |
| `nlm_k()` | Continuous XY point coordinates or CHILE point result | Translation-corrected K, centered L, annular g | `r`, `window` |
| `nlm_metrics()` | Categorical or thresholded raster | Per-patch and per-class metrics | `threshold`, `directions` |
| `nlm_compare()` | Observation and list of simulations | Observed curve, simulated median and pointwise quantiles | `statistic`, `r`, `angle`, `probs` |

### Directional spatial autocorrelation

```r
x <- nlm_anisotropic(ncol = 100, nrow = 100,
                     range_major = 35, range_minor = 6, angle = 45, seed = 1)
v <- nlm_variogram(x, lag = 3, max_dist = 45,
                   angle = c(0, 45, 90, 135), pairs = 100000, seed = 2)
plot(gamma ~ distance, v[v$angle == '45', ], type = 'b')
```

Each empirical variogram uses randomly sampled cell pairs. Values near large distances can have fewer sampled pairs. `angle` is measured counterclockwise from +x and is axial (0 and 180 are equivalent). Changing `seed` changes the sampled pairs.

### Point processes and observed stems

```r
obs <- nlm_points(ncol = 100, nrow = 100, method = 'thomas', n = 120, seed = 1)
k <- nlm_k(obs, r = seq(2, 25, by = 1))
plot(k$r, k$L, type = 'l', xlab = 'Distance', ylab = 'L(r) - r')

sims <- lapply(2:41, function(i) nlm_points(
  ncol = 100, nrow = 100, method = 'thomas', n = 120, seed = i))
comparison <- nlm_compare(obs, sims, statistic = 'L', r = seq(2, 25, by = 1))
head(comparison$curve)
```

For observed XY data, provide the same observation window used to generate the simulations, e.g. `window = c(0, 100, 0, 100)`. This initial `nlm_k()` handles **homogeneous** processes and rectangular observation windows only; spatially varying intensity is a potential confound. Use a point process package for inhomogeneous K and formal simultaneous global envelopes.

### Patch characteristics

```r
patches <- nlm_blob(ncol = 100, nrow = 100, n = 20, seed = 3)
summary <- nlm_metrics(patches)
summary$patches
summary$summary

# Thresholded spatial field (one foreground class)
x <- nlm_ambient(ncol = 100, nrow = 100, seed = 1)
nlm_metrics(x, threshold = 0.7)$summary
```

Patch metrics assume a planar equal-sized raster grid. `directions` determines connected components (4 or 8 neighbors); perimeter always counts exposed cell sides. Object IDs in a CHILE geometry result are treated as separate classes, so overlapping objects can fragment into multiple visible patches. For a single foreground class, pass a binary raster or threshold a continuous one.

**Interpretation:** With 20 simulations, pointwise 2.5% and 97.5% quantiles are descriptive and poorly resolved at the tails. Exceeding those curves is not automatically a formal model rejection. Comparing an observed pattern with fitted model simulations also requires accounting for parameter estimation and the observation-window design.
