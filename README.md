# CHILE

**Chaos in Silico Landscape Engine** · v0.9 development reference

CHILE generates synthetic landscapes from independent **random fields**, **point processes**, and **geometric objects**. You can simulate a landscape with no inputs, use an existing `terra::SpatRaster` as a template, or use an environmental raster to influence where objects or points are generated. Composition is a separate step.

CHILE is designed for controlled spatial experiments: vary autocorrelation, anisotropy, point clustering, patch geometry, or connectivity while holding other inputs fixed. The algorithms below describe the **current v0.9 source**, including its limitations—not hypothetical future features.

## Start here

```r
library(CHILE)
library(terra)   # attaches SpatRaster plot methods

nlm_list()       # available functions and families
nlm_show()       # default gallery

x <- nlm_anisotropic(range_major = 60, range_minor = 8, angle = 35, seed = 1)
plot(x)

p <- nlm_points(method = "thomas", n = 200, parents = 10,
                cluster_sd = 12, seed = 1)
plot(p$raster)
plot(p$points$x, p$points$y, asp = 1, pch = 16, cex = 0.5)

b <- nlm_blob(n = 40, environment = x, seed = 1)
plot(b$raster)
```

**Requirements:** R ≥ 4.1; `terra` is the spatial dependency. The package uses base R and the standard `stats`/`graphics` namespaces. No `sf`, DEM, D8, or external cost raster is needed. The standalone `CHILE` source package includes conventional `.Rd` help files, so `?nlm_points` should work in a correctly installed package. If R's Help pane returns an internal-server error, check `help("nlm_points", package="CHILE")` or `tools::Rd2txt(system.file("man/nlm_points.Rd", package="CHILE"))` as a diagnostic (the latter may depend on whether the installed package retains source `.Rd` files).

### Shared conventions

Most constructors accept `ncol = 256`, `nrow = 256`, `resolution = 1`, `template = NULL`, and `seed = NULL`. With no template, the output occupies `0..ncol*resolution` by `0..nrow*resolution`. With `template`, its grid and extent take precedence over the size arguments. The engine does not assign a real-world CRS by default. **Use projected map units** if distance, range, radius, or density are meant to represent meters or another linear unit.

`angle` is in degrees, nominally counterclockwise from the positive x-axis for field and point-process anisotropy. Raster geometry functions sometimes operate in row/column coordinates, so treat geometric orientation as visual until tested. `seed` controls a reproducible realization for most stochastic functions; `nlm_show()` defaults to `seed = 42`.

**Return types:** random-field functions return a `SpatRaster`, except `nlm_crosscorrelated()` which returns a two-layer `SpatRaster`. `nlm_points()` returns `list(raster, points, parameters)`, with **continuous subcell XY locations** and a cell-count raster. `nlm_polygon()`, `nlm_blob()`, and `nlm_bsp()` return `list(raster, objects)`; `nlm_filament()` returns `list(raster, paths)`; `nlm_partition()` returns `list(raster)`.

**Units matter:** covariance ranges are raster/map distances (with the FFT implementations interpreting them through cell resolution); point-process `intensity` is points per square map unit; `cluster_sd`, `radius`, `min_dist`, and `spacing` are map units. Polygon sizes and filament widths/lengths are **cells**, not meters. `n` means points or objects according to function; don't assume exact realized counts in every point process.

---

## 1. Random fields

### `nlm_random()` — independent uniform noise

```r
nlm_random(ncol=256, nrow=256, resolution=1,
           minval=0, maxval=1, template=NULL, seed=NULL)
```

- `minval`, `maxval`: lower and upper bounds of independent cell values.
- No autocorrelation or directional structure. Useful as null baseline and random perturbation.

```r
x <- nlm_random(minval = -1, maxval = 1, seed = 4)
```

### `nlm_gaussianfield()` — exponentially smooth Gaussian-style field

```r
nlm_gaussianfield(ncol=256, nrow=256, resolution=1,
                  range=20, sill=1, mean=0, angle=0,
                  range_minor=range, template=NULL, seed=NULL)
```

- `range`: major-axis spatial correlation scale.
- `range_minor`: minor-axis scale; set below `range` for elongated structure.
- `angle`: orientation of major axis.
- `mean`: output location; `sill`: output variance scale (`sd = sqrt(sill)`).
- Implemented through `nlm_anisotropic()` with Gaussian spectral filtering. Despite its historical name, this **is not an exponential-covariance implementation** in v0.9.

```r
x <- nlm_gaussianfield(range = 50, range_minor = 10,
                       angle = 45, sill = 2, seed = 1)
```

### `nlm_gaussianfield_matern()` — Matérn-like covariance

```r
nlm_gaussianfield_matern(ncol=256, nrow=256, resolution=1,
                         range=20, sill=1, mean=0, nu=1,
                         angle=0, range_minor=range,
                         template=NULL, seed=NULL)
```

- `range`, `range_minor`, `angle`: directional correlation scales and orientation.
- `nu`: smoothness; smaller values produce rougher fields, larger values smoother ones.
- `sill`, `mean`: marginal variance scale and mean.
- Uses a **spectral approximation** with periodic FFT boundaries, not exact Matérn simulation.

```r
rough  <- nlm_gaussianfield_matern(nu = 0.4, seed = 1)
smooth <- nlm_gaussianfield_matern(nu = 2.5, seed = 1)
nlm_show(rough = rough, smooth = smooth)
```

### `nlm_anisotropic()` — directional Gaussian field

```r
nlm_anisotropic(ncol=256, nrow=256, resolution=1,
                range_major=35, range_minor=8, angle=0,
                mean=0, sd=1, template=NULL, seed=NULL)
```

- `range_major`: correlation distance along the primary axis.
- `range_minor`: correlation distance across it.
- `angle`: orientation in degrees; 0 means along x, 90 along y.
- `mean`, `sd`: desired realized location and scale (standardized FFT realization).
- Ratio `range_major/range_minor` controls directional elongation. Equal ranges yield isotropic covariance. Ranges greater than the domain width can make realizations degenerate or dominated by broad structure.

```r
x <- nlm_anisotropic(range_major = 75, range_minor = 7,
                     angle = 30, seed = 10)
```

**Anisotropy in CHILE:** Gaussian and Matérn fields accept major/minor scales and angle; `nlm_mosaic()` accepts anisotropic *distance metrics*; `nlm_points(method="anisotropic_cluster")` uses oriented offspring clouds; `nlm_filament()` currently controls heading persistence rather than field-style covariance anisotropy. These do not all represent the same parameter.

### `nlm_multiscale()` — independent scales mixed together

```r
nlm_multiscale(ncol=256, nrow=256, resolution=1,
               ranges=c(8,30,90), weights=c(1,1,1),
               mean=0, sd=1, template=NULL, seed=NULL)
```

- `ranges`: vector of correlation scales.
- `weights`: signed contribution of each scale, same length as `ranges`; zero removes a component.
- `mean`, `sd`: output mean and standard deviation after recombination.
- Components are independent isotropic Gaussian fields; outputs are standardized after summing. Weight magnitudes set relative contribution, **not** exact realized variance fractions.

```r
x <- nlm_multiscale(ranges = c(4, 20, 80),
                    weights = c(0.3, 1, 2), seed = 3)
```

### `nlm_nonstationary()` — spatially changing texture

```r
nlm_nonstationary(ncol=256, nrow=256, resolution=1,
                  range_small=5, range_large=55,
                  transition=NULL, variance_gradient=0,
                  mean=0, sd=1, template=NULL, seed=NULL)
```

- `range_small`, `range_large`: scales of two independently generated component fields.
- `transition`: optional aligned raster of mixing weights (clamped to `[0,1]`); absent means a left-to-right transition.
- `variance_gradient`: `[-1,1]`, scales amplitude across the transition.
- `mean`, `sd`: overall baseline location and amplitude.
- The covariance changes across space; this is a two-scale spatial mixture, not a fully general nonstationary covariance model.

```r
x <- nlm_nonstationary(range_small = 3, range_large = 60,
                       variance_gradient = 0.5, seed = 5)
```

### `nlm_spectral()` — power-law spectral texture

```r
nlm_spectral(ncol=256, nrow=256, resolution=1,
             beta=2, low_cutoff=0.02,
             mean=0, sd=1, template=NULL, seed=NULL)
```

- `beta`: spectral slope. `0` approximates white noise; larger values concentrate energy at broad spatial frequencies.
- `low_cutoff`: regularizes the lowest frequencies; larger values suppress extremely broad variation.
- `mean`, `sd`: marginal scale.
- Stationary periodic Fourier synthesis; not a physical terrain model.

```r
nlm_show(fine = nlm_spectral(beta = 0, seed = 1),
         broad = nlm_spectral(beta = 4, seed = 1))
```

### `nlm_fbm()` — approximate fractional Brownian surface

```r
nlm_fbm(ncol=256, nrow=256, resolution=1,
        H=0.7, template=NULL, seed=NULL)
```

- `H`: Hurst-like roughness parameter, strictly between 0 and 1.
- Internally calls `nlm_spectral(beta = 2*H + 2)` and anchors the first cell at zero.
- **Approximation only**, not mathematically exact fractional Brownian motion with stationary increments.

### `nlm_midpoint()` — midpoint-like fractal texture

```r
nlm_midpoint(size=129, roughness=0.7,
             ncol=size, nrow=size, resolution=1,
             template=NULL, seed=NULL)
```

- `size`: convenient default for square grid dimensions; it is *not* restricted to `2^k + 1` in this release.
- `roughness`: strictly in `(0,1)`; controls weights in a four-scale spectral approximation.
- Despite the name, **this version does not implement the diamond-square midpoint-displacement algorithm**.

### `nlm_ambient()` — convenient noise presets

```r
nlm_ambient(ncol=256, nrow=256, resolution=1,
            extent=NULL, crs="", type="perlin",
            frequency=0.02, fractal=NULL, octaves=4,
            lacunarity=2, gain=0.5, normalize=TRUE,
            rescale_range=c(0,1), template=NULL, seed=NULL)
```

- `type`: **`"perlin"`, `"simplex"`, `"value"`, `"worley"`, `"white"`**.
- `frequency`: inverse characteristic scale; affects the range of component fields, or seed density in `worley`.
- `fractal`: `NULL`, `"billow"`, or `"rigidmulti"` have implemented transformations. Other strings are not validated and have no special effect.
- `octaves`: number of multiscale components.
- `lacunarity`: scale ratio across successive octaves.
- `gain`: successive octave weights.
- `normalize`: rescale the output to `rescale_range`; if `FALSE`, retain raw field values.
- `extent`, `crs`: optional spatial metadata when no `template` is provided.
- Important: `perlin`, `simplex`, and `value` are **currently the same underlying multiscale Gaussian approximation**. These are not distinct implementations of their named noise algorithms. `worley` is nearest-feature distance; `white` is independent uniform noise.

```r
x <- nlm_ambient(type = "worley", frequency = 0.015, seed = 2)
y <- nlm_ambient(type = "perlin", octaves = 5, gain = 0.4,
                 lacunarity = 2, seed = 2)
```

### `nlm_planargradient()` — linear environmental gradient

```r
nlm_planargradient(ncol=256, nrow=256, resolution=1,
                   angle=0, template=NULL)
```

- `angle`: gradient direction, degrees.
- Deterministic normalized output `[0,1]`; no `seed` because it uses no randomness.

### `nlm_radial()` — gradient around a center

```r
nlm_radial(ncol=256, nrow=256, resolution=1,
           center=NULL, decay=1, invert=FALSE, template=NULL)
```

- `center`: numeric XY center; default is center of raster extent.
- `decay`: positive exponent on normalized distance; below 1 spreads the gradient, above 1 concentrates change outward.
- `invert=TRUE`: largest values at the center instead of the edges.

### `nlm_distance()` — nearest-feature distance

```r
nlm_distance(ncol=256, nrow=256, resolution=1,
             n_features=10, features=NULL,
             invert=FALSE, template=NULL, seed=NULL)
```

- `n_features`: count of independently sampled reference points.
- `features`: optional two-column XY matrix replacing generated points.
- `invert`: reverse distance values as `max(distance) - distance`.
- Output values are distances in the raster's map units, **not** normalized `[0,1]`.

```r
x <- nlm_distance(n_features = 5, invert = TRUE, seed = 1)
```

### `nlm_percolation()` — independent binary occupancy

```r
nlm_percolation(ncol=256, nrow=256, resolution=1,
                prob=0.5, template=NULL, seed=NULL)
```

- `prob`: probability each cell independently equals 1, from 0 to 1.
- Despite the name, this is an independent Bernoulli lattice; it does not force a connected or spanning cluster.

### `nlm_neighborhood()` — correlated binary occupancy

```r
nlm_neighborhood(ncol=256, nrow=256, resolution=1,
                 beta=0.7, intercept=0, sweeps=25,
                 neighbors=4, initial_prob=0.5,
                 template=NULL, seed=NULL)
```

- `beta`: neighbor coupling strength. Positive values promote agreement/patchiness; negative values favor alternation.
- `intercept`: baseline bias toward occupied cells. It is **not** an exact target occupancy proportion.
- `sweeps`: number of full random-order Gibbs updates; too few may retain initialization effects.
- `neighbors`: **4 or 8** cell neighbors.
- `initial_prob`: initial independent occupancy probability.
- This is an Ising-like Gibbs sampler, not a controlled-patch-size generator.

```r
x <- nlm_neighborhood(beta = 1.1, intercept = -0.3,
                      sweeps = 60, neighbors = 8, seed = 1)
```

### `nlm_curdle()` — hierarchical binary subdivisions

```r
nlm_curdle(ncol=256, nrow=256, resolution=1,
           levels=6, p=0.7, template=NULL, seed=NULL)
```

- `levels`: number of recursive quartering steps; `0` returns all ones.
- `p`: probability each child of an active rectangle remains active.
- Large `levels` gives finer hierarchical detail; small `p` produces sparse survivors.
- Unlike `nlm_partition()`, this creates binary occupancy, not polygon IDs.

### `nlm_crosscorrelated()` — paired continuous fields

```r
nlm_crosscorrelated(ncol=256, nrow=256, resolution=1,
                    rho=0.6, range1=30, range2=10,
                    template=NULL, seed=NULL)
```

- `rho`: intended **realized zero-lag Pearson correlation** between two layers, from -1 to 1.
- `range1`, `range2`: the first field's scale and the second independent component's scale.
- Returns a **two-layer** raster: `field1`, `field2`.
- Orthogonalization enforces the intended correlation in the finite sample; the second field inherits structure from the first. This is not a general multivariate covariance model.

```r
xy <- nlm_crosscorrelated(rho = -0.7, range1 = 60,
                          range2 = 8, seed = 3)
plot(xy)
```

### `nlm_mosaic()` — anisotropic Voronoi mosaic

```r
nlm_mosaic(ncol=256, nrow=256, resolution=1,
           npts=20, angle=0, range_major=1,
           range_minor=1, template=NULL, seed=NULL)
```

- `npts`: seed count, creating nearest-seed regions.
- `angle`: rotation of anisotropic distance metric.
- `range_major`, `range_minor`: directional scaling of distance to seeds. Their **ratio** changes cell shapes; multiplying both by the same factor does not.
- Returns a categorical raster of region IDs, not point coordinates.

```r
x <- nlm_mosaic(npts = 35, range_major = 4,
                range_minor = 1, angle = 45, seed = 1)
```

`nlm_mosaiccluster(npts=20, cluster_sd=5, ...)` is retained as an alias to `nlm_mosaic()`. In this version **`cluster_sd` does nothing**; use `nlm_mosaic()` directly.

---

## 2. Point processes: `nlm_points()`

This is the primary interface for simulating individual locations, events, stems, or object origins. Coordinates are stored as continuous XY rather than rounded to cells. A companion raster counts how many locations fall in each cell.

```r
nlm_points(ncol=256, nrow=256, resolution=1,
  method="binomial", n=100, intensity=NULL, environment=NULL,
  parents=12, offspring=10, cluster_sd=12, radius=15,
  min_dist=8, interaction=0.5,
  range_major=35, range_minor=8, angle=0,
  spacing=NULL, jitter=0.3,
  mixture=c("thomas","inhibition"), weights=NULL,
  max_attempts=100000, template=NULL, seed=NULL)
```

### Available methods

| `method` | Model / behavior | Main controls | Important qualification |
|---|---|---|---|
| `"binomial"` | Fixed-count uniform random points | `n` | Ordinary CSR conditioned on count |
| `"poisson"` | Homogeneous Poisson locations | `n=NULL`, `intensity` | With `n` supplied, becomes count-conditioned (binomial) |
| `"inhomogeneous"` | Independently located points weighted by a surface | `environment`, `n` or `intensity` | If no `environment`, behaves like homogeneous placement |
| `"thomas"` | Gaussian offspring around random parents | `parents`, `offspring`, `cluster_sd`, `n` | Parent count is fixed unless `n=NULL` |
| `"matern_cluster"` | Offspring uniformly in parent-centered discs | `parents`, `offspring`, `radius`, `n` | `radius` is cluster radius |
| `"anisotropic_cluster"` | Elliptical Gaussian parent-offspring clouds | `parents`, `range_major`, `range_minor`, `angle`, `n` | Scales are offspring standard deviations |
| `"inhibition"` | Sequentially accept points above minimum distance | `min_dist`, `max_attempts`, `n` | Can fail to reach requested `n` |
| `"strauss"` | Approximately repulsive Gibbs interaction | `radius`, `interaction`, `max_attempts`, `n` | Fixed-N approximate Metropolis sampler |
| `"lgcp"` | Poisson / count-conditioned points from latent log-Gaussian intensity | `range_major`, `range_minor`, `angle`, `environment`, `n` | Simplified latent field; intensity is normalized to target expected count |
| `"grid"` | Regularly spaced points | `spacing`, `angle`, `n` | May produce fewer than requested points |
| `"jittered"` | Randomly displaced regular grid | `spacing`, `jitter`, `angle`, `n` | `jitter` from 0 to 1 |
| `"mixture"` | Combined independently generated components | `mixture`, `weights`, `n` | Shares parameter values across component models |

### Point parameters, explained

**Count and density**

- `n=100`: requested point count. For binomial/Poisson-conditioned processes, it is the realized count. For inhibition, grids, and some cluster realizations it may not be achieved.
- `n=NULL`: use process-derived stochastic counts. For Poisson-like cases, count is drawn from a Poisson law with expected value `intensity * area`; if `intensity=NULL`, the code uses an expected count of 100. For cluster models the parent and offspring counts are stochastic.
- `intensity=NULL`: point density **per square map unit**, not per raster cell. With a large projected extent, `intensity=0.1` can imply enormous counts.

**Environmental influence**

- `environment=NULL`: uniform origin placement or internally simulated latent fields as appropriate.
- `environment=<SpatRaster>`: aligned **nonnegative** weights used for origin locations / relative intensity for supported modes. High cell values attract more origins; the surface does **not** change subsequent offspring shapes or repulsive interactions. Grid modes do not use environmental weighting. For `lgcp`, the environment multiplies the latent intensity.
- `template=NULL`: sets grid only, not placement weights. When provided, it must be aligned with `environment`.

**Clustering**

- `parents=12`: number of cluster centers when count is fixed; for `n=NULL`, the code draws a Poisson number of parents with this expectation.
- `offspring=10`: mean offspring **per parent**, used when `n=NULL`; when `n` is supplied, the requested number is allocated among parents instead.
- `cluster_sd=12`: standard deviation of Gaussian offsets for `thomas`.
- `radius=15`: uniform-disc radius for `matern_cluster`; Strauss interaction radius for `strauss`.
- `range_major=35`, `range_minor=8`, `angle=0`: covariance ellipse for `anisotropic_cluster`, and anisotropic latent intensity texture in `lgcp`.

**Repulsion / organization**

- `min_dist=8`: minimum interpoint distance for sequential `inhibition`. Increase it to create increasingly regular patterns; ensure the domain can accommodate the requested count.
- `interaction=0.5`: Strauss pair interaction parameter, constrained to `[0,1]`; `1` means no pair penalty, `0` strongly rejects configurations that create additional within-radius pairs. This implementation does **not** support an attractive Strauss model (`interaction > 1`).
- `spacing=NULL`: approximate lattice spacing computed from area and requested count, or user-specified in map units.
- `jitter=0.3`: fraction of lattice spacing for independent XY perturbations; only used in `jittered`.
- `max_attempts=100000`: stopping limit for rejection/inhibition and cluster edge correction; Metropolis updates in Strauss are capped at `min(max_attempts, 50*n)`.

**Mixtures**

- `mixture=c("thomas","inhibition")`: component method names, excluding `"mixture"` itself.
- `weights=NULL`: equal proportions of the requested count among components, or provide nonnegative weights. These are **allocation weights**, not mixture-field covariance weights.
- All components receive the same values for common parameters such as `radius`, `min_dist`, and `angle`. Component-specific settings are **not yet supported**.

### Point-process recipes

```r
# Exactly 500 independently located points
p1 <- nlm_points(method = "binomial", n = 500, seed = 1)

# True random Poisson count: expected 500 points on a 200 x 200 domain
p2 <- nlm_points(ncol = 200, nrow = 200, method = "poisson",
                 n = NULL, intensity = 500 / (200 * 200), seed = 1)

# Tight Gaussian groups
p3 <- nlm_points(method = "thomas", n = 500,
                 parents = 20, cluster_sd = 5, seed = 1)

# Broad circular colonies
p4 <- nlm_points(method = "matern_cluster", n = 500,
                 parents = 12, radius = 25, seed = 1)

# Elongated groups
p5 <- nlm_points(method = "anisotropic_cluster", n = 500,
                 parents = 20, range_major = 35,
                 range_minor = 4, angle = 45, seed = 1)

# Spaced-out stems / trees
p6 <- nlm_points(method = "inhibition", n = 200,
                 min_dist = 10, seed = 1)

# Soft rather than absolute repulsion
p7 <- nlm_points(method = "strauss", n = 200,
                 radius = 12, interaction = 0.2, seed = 1)

# Grouping driven by an internally simulated intensity field
p8 <- nlm_points(method = "lgcp", n = 400,
                 range_major = 50, range_minor = 10,
                 angle = 30, seed = 1)

# Planted grid with position jitter
p9 <- nlm_points(method = "jittered", n = 250,
                 spacing = 14, jitter = 0.5,
                 angle = 20, seed = 1)

# Independent environmental field controls origin probability
env <- nlm_anisotropic(ncol = 200, nrow = 200, seed = 1)
w <- terra::app(env, function(z) exp(z))
p10 <- nlm_points(method = "inhomogeneous", n = 400,
                  environment = w, seed = 1)

# Mix clustered and inhibited patterns
p11 <- nlm_points(method = "mixture", n = 300,
                  mixture = c("thomas", "inhibition"),
                  weights = c(0.7, 0.3), seed = 1)

nlm_show(CSR = p1, clustered = p3, inhibited = p6, organized = p9)
```

**Interpreting point output:** `p$points` stores `id`, `x`, `y`; `p$raster` stores per-cell **counts**, not binary presence; `p$parameters` stores a short subset of requested settings and realized count. Edge clipping and sequential inhibition can reduce count. `n` for grid/jittered is a cap or spacing-derived target, not a guarantee. Point-process realizations have boundary effects; no toroidal edge correction is used.

---

## 3. Spatial geometry

Geometry functions operate without a field. `environment` is an optional **placement-weighting surface**: it selects object centers, not shapes, unless noted. All geometry rasters use object IDs, with `NA` outside generated objects.

### `nlm_polygon()` — straight-edged patches

```r
nlm_polygon(ncol=256, nrow=256, resolution=1,
            n=30, min_size=8, max_size=50, max_ratio=4,
            angle=c(0,180), skew=0, taper=0,
            environment=NULL, template=NULL,
            overwrite=TRUE, seed=NULL)
```

- `n`: number of proposed polygons.
- `min_size`, `max_size`: range of sampled widths (in cells).
- `max_ratio`: width-to-height ratio drawn uniformly between 1 and this maximum.
- `angle`: **two-number interval** of rotation angles in degrees, e.g. `c(20,60)`; use `c(45,45)` for fixed rotation.
- `skew`: top-edge horizontal offset, in units proportional to width.
- `taper`: scales the top two vertices; `0` maintains full width. Large values may create self-intersections.
- `environment`: weights center placement.
- `overwrite`: newer polygons replace earlier IDs where they overlap; `FALSE` keeps earlier object IDs.
- Returns `list(raster, objects)` with center and sampled shape parameters. Actual visible areas differ from nominal areas because of clipping and overlap.

```r
p <- nlm_polygon(n = 60, min_size = 20, max_size = 70,
                 angle = c(20, 65), skew = 0.3,
                 taper = 0.25, seed = 1)
plot(p$raster)
```

### `nlm_blob()` — irregular elliptical patches

```r
nlm_blob(ncol=256, nrow=256, resolution=1,
         n=30, min_size=8, max_size=50, max_ratio=4,
         irregularity=0.3, complexity=5,
         environment=NULL, template=NULL,
         overwrite=TRUE, seed=NULL)
```

- `n`: proposed blobs.
- `min_size`, `max_size`: sampled major-axis diameters in cells.
- `max_ratio`: maximum elongation of ellipse.
- `irregularity`: amplitude of sinusoidal perturbation to elliptical boundary; `0` yields ellipses. High values can make shapes highly distorted.
- `complexity`: number of angular harmonic terms in the boundary perturbation; higher values give finer scalloping.
- `environment`, `overwrite`: same meaning as in polygons.
- Returns `list(raster, objects)`. Boundary perturbation changes area; prescribed area is **not** enforced.

```r
b <- nlm_blob(n = 50, irregularity = 0.65,
              complexity = 9, max_ratio = 5, seed = 3)
```

### `nlm_filament()` — linear corridors and networks

```r
nlm_filament(ncol=256, nrow=256, resolution=1,
             n=20, length=c(30,100), width=c(1,5),
             mode="walk", tortuosity=0.2, branching=0.02,
             environment=NULL, template=NULL, seed=NULL)
```

`mode` choices:

| Mode | Process | Controls / caution |
|---|---|---|
| `"straight"` | Persistent random heading | `length`, `width`; straight lines |
| `"walk"` | Correlated random heading | `tortuosity` is SD of heading increments in **radians per step** |
| `"branching"` | Correlated walk with secondary branches | `tortuosity`, `branching` is **per-step** fork probability; capped at 100 segments per origin |
| `"field"` | Step toward the lowest-valued neighboring cell | `environment` supplies elevation/cost-like raster; if `NULL`, internally generates a multiscale field |

- `n`: number of initial origins, not the final number of branches.
- `length=c(min,max)`: sampled maximum number of path steps in cells, not guaranteed final length (edges can truncate).
- `width=c(min,max)`: sampled raster corridor diameter, in cells.
- `tortuosity`: zero yields nearly straight movement in walk/branching; increasing it increases heading variability.
- `branching`: used only when `mode="branching"`; probability evaluated at each step.
- `environment`: weights starting cell selection in all modes. In `"field"` it also guides steps toward local minima.
- Returns `list(raster, paths)`; `paths` contains matrices of **row/column cell coordinates**, grouped by initial origin and segment.

**Current limitation:** `"field"` is a local minimum-following heuristic, **not** a least-cost path algorithm, true hydrological drainage routing, or continuous-gradient tracing. It can stop in a local minimum. `nlm_filament()` currently does not expose `overwrite`; later objects overwrite earlier IDs.

```r
f1 <- nlm_filament(mode = "straight", n = 15, seed = 1)
f2 <- nlm_filament(mode = "walk", n = 15, tortuosity = 0.5, seed = 1)
f3 <- nlm_filament(mode = "branching", n = 15,
                   tortuosity = 0.1, branching = 0.04, seed = 1)
f4 <- nlm_filament(mode = "field", n = 15, seed = 1)
nlm_show(straight = f1, winding = f2, branches = f3, guided = f4)
```

### `nlm_partition()` — complete landscape subdivisions

```r
nlm_partition(ncol=256, nrow=256, resolution=1,
              n=25, method="voronoi", angle=0,
              template=NULL, seed=NULL)
```

- `method="voronoi"`: nearest-seed regions, uses `nlm_mosaic(npts=n)`.
- `method="grid"`: rotated square lattice; `angle` applies here, and `n` is an approximate target (not exact number of cells/polygons).
- `method="bsp"`: recursively split the largest rectangle into two axis-aligned rectangles until reaching `n` or the smallest possible cells. This is the **true recursive BSP** implementation, unlike `nlm_bsp()`.
- Returns `list(raster)` with IDs occupying the full grid.

### `nlm_bsp()` — convenience alias for rectangles

```r
nlm_bsp(ncol=256, nrow=256, resolution=1,
        n_patches=100, min_size=10, max_size=80,
        max_ratio=4, environment=NULL, template=NULL,
        overwrite=TRUE, seed=NULL)
```

Calls `nlm_polygon()` with `angle=c(0,0)`. This is a **rectangle-stamping** generator, not recursive binary space partitioning. For the latter use `nlm_partition(method="bsp")`.

---

## 4. Transform, combine, compose

### `nlm_transform()` — choose marginal distribution

```r
nlm_transform(nlm,
  distribution=c("normal","uniform","lognormal",
                 "gamma","beta","exponential"),
  params=NULL, zero_prob=0, seed=NULL)
```

- `distribution`: **`"normal"`, `"uniform"`, `"lognormal"`, `"gamma"`, `"beta"`, `"exponential"`**.
- `params`: named list passed to the corresponding R quantile function. Defaults: normal `mean=0,sd=1`; uniform `min=0,max=1`; lognormal `meanlog=0,sdlog=1`; gamma `shape=2,scale=1`; beta `shape1=2,shape2=2`; exponential `rate=1`.
- `zero_prob`: fraction of finite cells randomly reset to exactly zero (ties introduced).
- Rank-based quantile mapping preserves **rank ordering**, not Pearson covariance or exact variogram values.

```r
x <- nlm_anisotropic(seed = 1)
y <- nlm_transform(x, "gamma", params = list(shape = 2, scale = 3))
z <- nlm_transform(x, "lognormal", zero_prob = 0.3, seed = 1)
```

### `nlm_mix()` — weighted mixture of continuous fields

```r
nlm_mix(..., weights=NULL, mean=0, sd=1)
```

- `...`: one or more aligned, complete `SpatRaster` fields.
- `weights`: one weight per field; defaults equal. Each field is individually standardized before combining.
- `mean`, `sd`: rescale the final combined field; weights are relative contributions, not guaranteed output variance fractions.
- Requires complete finite values. A single constant field can yield nonfinite values when standardized.

```r
x <- nlm_anisotropic(seed = 1)
y <- nlm_spectral(seed = 2)
z <- nlm_mix(x, y, weights = c(0.8, 0.2))
```

### `nlm_compose()` — overlay discrete objects on a landscape

```r
nlm_compose(background, objects,
            method=c("replace","add","multiply"),
            transparent=NA_real_)
```

- `background`: continuous/categorical raster used as base.
- `objects`: aligned `SpatRaster` or geometry result containing `$raster`.
- `method`: **`"replace"`**, **`"add"`**, **`"multiply"`**.
- `transparent`: object value to exclude, with `NA` meaning all non-finite object cells are ignored.
- Object **ID values** are used directly in arithmetic operations. If you mean a binary footprint, convert to 0/1 first; otherwise `multiply` multiplies by object IDs, which is rarely the desired ecological operation.

```r
bg <- nlm_multiscale(ncol = 200, nrow = 200, seed = 1)
ob <- nlm_blob(n = 30, template = bg, seed = 1)
composite <- nlm_compose(bg, ob, method = "replace")

# Binary treatment-mask example
mask <- !is.na(ob$raster)
treated <- terra::ifel(mask, 0, bg)  # set treated cells to zero
```

---

## 5. Discover and visualize: `nlm_list()` and `nlm_show()`

```r
nlm_list(family=NULL)
```

- `family`: optional `"field"`, `"points"`, `"geometry"`, `"operation"`, or `"utility"`.
- Returns a data frame of documented function names and families.

```r
nlm_show(..., methods=NULL, family=NULL,
         ncol=90, nrow=90, seed=42)
```

- With **no arguments**, generates a 12-panel gallery of representative generators.
- `methods`: character vector of function names. Each function is called with its defaults on `ncol` × `nrow` grid, except that the gallery uses `nlm_points(method="thomas", n=150)` and `nlm_filament(n=12)`.
- `family`: optional filter on method names when `methods` is supplied.
- `...`: any named `SpatRaster` or a result list containing `$raster`, e.g. point or geometry output. Best way to compare custom parameters.
- `seed`: deterministic gallery seed; restores the prior RNG state after running methods.
- Returns the displayed results as an **invisible list**.

```r
nlm_show()
nlm_show(methods = c("nlm_anisotropic", "nlm_blob", "nlm_filament"))

x <- nlm_anisotropic(range_major = 60, range_minor = 5, angle = 0, seed = 1)
y <- nlm_anisotropic(range_major = 60, range_minor = 5, angle = 90, seed = 1)
nlm_show(horizontal = x, vertical = y)

nlm_list("points")
```

If `plot(x)` unexpectedly dispatches to `plot.default()` after `library(CHILE)`, attach `library(terra)` or call `terra::plot(x)`; `CHILE` currently imports the `terra` namespace without attaching all of its methods.

---

## 6. Design patterns for experiments

### Separate pattern from object placement

```r
# One synthetic environment, three independently drawn feature geometries
env <- nlm_anisotropic(ncol = 250, nrow = 250,
                       range_major = 50, range_minor = 10,
                       angle = 35, seed = 1)

# Shift field values to positive relative placement weights
w <- env - terra::global(env, "min", na.rm = TRUE)[1,1] + 0.01

rectangles <- nlm_polygon(n = 40, environment = w, seed = 2)
blobs <- nlm_blob(n = 40, environment = w, seed = 2)
filaments <- nlm_filament(n = 40, environment = w,
                          mode = "walk", seed = 2)

nlm_show(rectangles = rectangles, blobs = blobs, filaments = filaments)
```

### Test anisotropy without changing other model settings

```r
angles <- c(0, 30, 60, 90)
fields <- lapply(angles, function(a)
  nlm_anisotropic(range_major = 50, range_minor = 8,
                  angle = a, seed = 2))
names(fields) <- paste0(angles, " degrees")
do.call(nlm_show, fields)
```

### Sample a continuous field at simulated point locations

```r
x <- nlm_multiscale(ncol = 200, nrow = 200, seed = 1)
p <- nlm_points(method = "thomas", n = 500, template = x, seed = 2)
# terra::extract takes cell numbers or a spatial vector; use cellFromXY here
cells <- terra::cellFromXY(x, as.matrix(p$points[, c("x", "y")]))
p$points$environment <- terra::values(x)[cells, 1]
```

---

## 7. Current development limitations

CHILE v0.9 is an **experimental implementation**, not a validated stochastic simulation library. Several labels represent useful approximations rather than exact named processes. In particular:

1. `nlm_gaussianfield()` uses Gaussian spectral smoothing rather than exponential covariance; Matérn is an FFT spectral approximation.
2. `nlm_fbm()` and `nlm_midpoint()` are spectral surrogates, not exact classical algorithms.
3. `nlm_ambient(type="perlin"/"simplex"/"value")` currently routes to the same Gaussian multiscale implementation.
4. `nlm_points(method="strauss")` uses finite-run Metropolis relocation, not perfect simulation; `"inhibition"` is sequential acceptance.
5. `nlm_filament(mode="field")` is local-neighbor descent, not validated drainage or optimal routing.
6. Polygon and blob dimensions are measured in cells; their area and visible count can change after clipping and overlap.
7. FFT fields have periodic boundary assumptions; finite realized spatial statistics can deviate from model parameters.
8. The source has not been executed through a full R validation matrix in this development environment. Inspect and test results before using them in inference.

**Recommended validation for v1.0:** realized marginal distributions; directional variograms and anisotropy ratios; spectral slopes; cluster pair-correlation and nearest-neighbor distances; point-count expectations; object areas/edge densities; reproducibility; and performance scaling with raster size and number of simulated objects.

---

## Function index

**Continuous or binary fields:** `nlm_random`, `nlm_gaussianfield`, `nlm_gaussianfield_matern`, `nlm_anisotropic`, `nlm_multiscale`, `nlm_nonstationary`, `nlm_spectral`, `nlm_fbm`, `nlm_midpoint`, `nlm_ambient`, `nlm_planargradient`, `nlm_radial`, `nlm_distance`, `nlm_percolation`, `nlm_neighborhood`, `nlm_curdle`, `nlm_crosscorrelated`, `nlm_mosaic`, `nlm_mosaiccluster`.

**Points:** `nlm_points` (12 methods).

**Geometry:** `nlm_polygon`, `nlm_blob`, `nlm_filament`, `nlm_partition`, `nlm_bsp`.

**Operations and utilities:** `nlm_transform`, `nlm_mix`, `nlm_compose`, `nlm_show`, `nlm_list`.
