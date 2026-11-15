# CHILE v0.9 — additional standalone neutral landscape generators
# Required: terra. Base R otherwise. No sf, sp, gstat, fields, or igraph.
# All dimensions and ranges are in raster cells unless a template is provided.
# In a supplied template, ranges are measured in map units (assuming projected CRS).
# Every generator runs without an input environmental raster.

# --------------------------- internal helpers ----------------------------
.chile_raster <- function(ncol = 256, nrow = 256, resolution = 1,
                          template = NULL) {
  if (!is.null(template)) {
    if (!inherits(template, "SpatRaster")) stop("template must be a SpatRaster")
    r <- template[[1]]
    terra::values(r) <- NA_real_
  } else {
    r <- terra::rast(ncols = ncol, nrows = nrow,
                     xmin = 0, xmax = ncol * resolution,
                     ymin = 0, ymax = nrow * resolution)
  }
  r
}

.chile_out <- function(mat, r, name) {
  terra::values(r) <- as.vector(t(mat))
  names(r) <- name
  r
}

.chile_standardize <- function(z) {
  z <- z - mean(z)
  s <- stats::sd(as.vector(z))
  if (is.finite(s) && s > 0) z <- z / s
  z
}

.chile_freq <- function(n) {
  k <- 0:(n - 1L)
  ifelse(k <= floor(n / 2), k, k - n) / n
}

# Generate a real stationary Gaussian field by filtering white noise in
# Fourier space. Spectral amplitudes are nonnegative and symmetric, so the
# inverse transform is real. Periodic boundary is an FFT approximation.
.chile_fft_field <- function(nr, nc, amplitude) {
  if (!is.matrix(amplitude) || !all(dim(amplitude) == c(nr, nc)))
    stop("amplitude must have dimensions nrow x ncol")
  w <- matrix(stats::rnorm(nr * nc), nr, nc)
  z <- Re(fft(fft(w) * amplitude, inverse = TRUE)) / (nr * nc)
  .chile_standardize(z)
}

# ---------------------- 1. anisotropic Gaussian --------------------------
# range_major / range_minor: Gaussian correlation length in raster units.
# angle is counterclockwise from positive x (map coordinates).
nlm_anisotropic <- function(ncol = 256, nrow = 256, resolution = 1,
                            range_major = 35, range_minor = 8,
                            angle = 0, mean = 0, sd = 1,
                            template = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (range_major <= 0 || range_minor <= 0 || sd < 0)
    stop("ranges must be positive and sd nonnegative")
  r <- .chile_raster(ncol, nrow, resolution, template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  fx <- matrix(rep(.chile_freq(nc) / terra::res(r)[1], each = nr), nr, nc)
  fy <- matrix(rep(.chile_freq(nr) / terra::res(r)[2], times = nc), nr, nc)
  a <- -angle * pi / 180
  u <- fx * cos(a) + fy * sin(a)
  v <- -fx * sin(a) + fy * cos(a)
  amp <- exp(-pi^2 * ((range_major * u)^2 + (range_minor * v)^2))
  z <- .chile_fft_field(nr, nc, amp)
  .chile_out(mean + sd * z, r, "anisotropic")
}

# ----------------------- 2. multiscale synthesis -------------------------
# Independent Gaussian fields mixed with user-specified scale weights.
nlm_multiscale <- function(ncol = 256, nrow = 256, resolution = 1,
                           ranges = c(8, 30, 90), weights = c(1, 1, 1),
                           mean = 0, sd = 1, template = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (length(ranges) != length(weights) || !length(ranges) ||
      any(ranges <= 0) || any(!is.finite(weights)) ||
      all(weights == 0) || sd < 0) stop("Invalid ranges, weights, or sd")
  r <- .chile_raster(ncol, nrow, resolution, template)
  z <- matrix(0, terra::nrow(r), terra::ncol(r))
  for (i in seq_along(ranges)) {
    if (weights[i] == 0) next
    f <- nlm_anisotropic(range_major = ranges[i], range_minor = ranges[i],
                         template = r)
    z <- z + weights[i] * matrix(terra::values(f)[, 1], nrow(z), ncol(z), byrow = TRUE)
  }
  .chile_out(mean + sd * .chile_standardize(z), r, "multiscale")
}

# ------------------------ 4. nonstationary field -------------------------
# Mixture of two spatial scales with a spatially varying transition.
# Gradient may be supplied by the user (0..1) or generated left-to-right.
# Both scales remain present near the midpoint; sd controls local amplitude.
nlm_nonstationary <- function(ncol = 256, nrow = 256, resolution = 1,
                              range_small = 5, range_large = 55,
                              transition = NULL, variance_gradient = 0,
                              mean = 0, sd = 1,
                              template = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (range_small <= 0 || range_large <= 0 || sd < 0 ||
      abs(variance_gradient) > 1) stop("Invalid parameters")
  r <- .chile_raster(ncol, nrow, resolution, template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  if (is.null(transition)) {
    t <- matrix(rep(seq(0, 1, length.out = nc), each = nr), nr, nc)
  } else {
    if (!terra::compareGeom(r, transition, stopOnError = FALSE))
      stop("transition must match the output grid")
    t <- matrix(terra::values(transition)[, 1], nr, nc, byrow = TRUE)
    t <- pmin(1, pmax(0, t))
  }
  a <- nlm_anisotropic(range_major = range_small, range_minor = range_small,
                       template = r)
  b <- nlm_anisotropic(range_major = range_large, range_minor = range_large,
                       template = r)
  za <- matrix(terra::values(a)[, 1], nr, nc, byrow = TRUE)
  zb <- matrix(terra::values(b)[, 1], nr, nc, byrow = TRUE)
  z <- (1 - t) * za + t * zb
  # Correct local mixture variance under independent standardized fields.
  z <- z / sqrt((1 - t)^2 + t^2)
  z <- z * (1 + variance_gradient * (2 * t - 1))
  .chile_out(mean + sd * z, r, "nonstationary")
}

# ------------------------ 5. neighborhood / Ising ------------------------
# Binary Gibbs random field. beta controls clustering, intercept controls
# prevalence. Sweeps perform randomly ordered single-cell Gibbs updates.
nlm_neighborhood <- function(ncol = 256, nrow = 256, resolution = 1,
                             beta = 0.7, intercept = 0, sweeps = 25,
                             neighbors = 4, initial_prob = 0.5,
                             template = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (!neighbors %in% c(4, 8) || sweeps < 0 ||
      initial_prob < 0 || initial_prob > 1) stop("Invalid parameters")
  r <- .chile_raster(ncol, nrow, resolution, template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  z <- matrix(stats::rbinom(nr * nc, 1, initial_prob), nr, nc)
  offsets <- rbind(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
  if (neighbors == 8)
    offsets <- rbind(offsets, c(-1,-1), c(-1,1), c(1,-1), c(1,1))
  for (s in seq_len(sweeps)) {
    for (cell in sample.int(nr * nc)) {
      row <- (cell - 1L) %% nr + 1L
      col <- (cell - 1L) %/% nr + 1L
      rr <- row + offsets[, 1]; cc <- col + offsets[, 2]
      ok <- rr >= 1 & rr <= nr & cc >= 1 & cc <= nc
      signal <- sum(2 * z[cbind(rr[ok], cc[ok])] - 1)
      z[row, col] <- stats::rbinom(1, 1,
        stats::plogis(intercept + beta * signal))
    }
  }
  .chile_out(z, r, "neighborhood")
}

# ------------------------- 6. hierarchical curdling ----------------------
# Each retained rectangle splits into four children, each kept with p.
# levels=1 is coarse; larger values yield fine-grained hierarchy.
nlm_curdle <- function(ncol = 256, nrow = 256, resolution = 1,
                       levels = 6, p = 0.7, template = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (levels < 0 || p < 0 || p > 1) stop("Invalid levels or p")
  r <- .chile_raster(ncol, nrow, resolution, template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  z <- matrix(0L, nr, nc)
  active <- list(c(1L, nr, 1L, nc))
  for (level in seq_len(levels)) {
    next_active <- list()
    for (rect in active) {
      r1 <- rect[1]; r2 <- rect[2]; c1 <- rect[3]; c2 <- rect[4]
      rm <- floor((r1 + r2) / 2); cm <- floor((c1 + c2) / 2)
      children <- list(
        c(r1, rm, c1, cm), c(r1, rm, cm + 1L, c2),
        c(rm + 1L, r2, c1, cm), c(rm + 1L, r2, cm + 1L, c2)
      )
      for (child in children) {
        if (child[1] > child[2] || child[3] > child[4]) next
        if (stats::runif(1) < p) {
          z[child[1]:child[2], child[3]:child[4]] <- 1L
          next_active[[length(next_active) + 1L]] <- child
        }
      }
    }
    active <- next_active
    if (!length(active)) break
  }
  # For levels=0, the starting region remains occupied.
  if (levels == 0) z[,] <- 1L
  .chile_out(z, r, "curdle")
}

# ---------------------- 7. correlated Gaussian fields -------------------
# rho is zero-lag Pearson correlation (approximately, in finite samples).
# Second field inherits part of first field's spatial scale. This is a
# valid construction; arbitrary separate covariance/cross-covariance
# specifications are not in general jointly admissible.
nlm_crosscorrelated <- function(ncol = 256, nrow = 256, resolution = 1,
                                rho = 0.6, range1 = 30, range2 = 10,
                                template = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (abs(rho) > 1) stop("rho must be between -1 and 1")
  r <- .chile_raster(ncol, nrow, resolution, template)
  a <- nlm_anisotropic(range_major = range1, range_minor = range1, template = r)
  b <- nlm_anisotropic(range_major = range2, range_minor = range2, template = r)
  x <- terra::values(a)[, 1]
  y <- terra::values(b)[, 1]
  x <- as.vector(scale(x))
  # Orthogonalize to deliver the requested finite-sample correlation.
  y <- y - mean(y)
  y <- y - sum(x * y) / sum(x^2) * x
  sy <- stats::sd(y)
  if (!is.finite(sy) || sy == 0) stop("Degenerate independent field")
  y <- y / sy
  z <- rho * x + sqrt(1 - rho^2) * y
  r <- c(r, r)
  terra::values(r) <- cbind(x, z)
  names(r) <- c("field1", "field2")
  r
}

# ----------------------- 8. spectral random field ------------------------
# Power spectrum proportional to (k^2 + k0^2)^(-beta/2).
# beta=0 gives white noise; higher beta favors broad features.
# Spectral amplitude is square root of the power spectrum.
nlm_spectral <- function(ncol = 256, nrow = 256, resolution = 1,
                         beta = 2, low_cutoff = 0.02,
                         mean = 0, sd = 1,
                         template = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (beta < 0 || low_cutoff <= 0 || sd < 0) stop("Invalid parameters")
  r <- .chile_raster(ncol, nrow, resolution, template)
  nr <- terra::nrow(r); nc <- terra::ncol(r)
  fx <- matrix(rep(.chile_freq(nc), each = nr), nr, nc)
  fy <- matrix(rep(.chile_freq(nr), times = nc), nr, nc)
  amp <- (fx^2 + fy^2 + low_cutoff^2)^(-beta / 4)
  z <- .chile_fft_field(nr, nc, amp)
  .chile_out(mean + sd * z, r, "spectral")
}

# ------------------------ 9. radial gradient -----------------------------
# Distance to a point in the raster's map units; normalized to 0..1.
nlm_radial <- function(ncol = 256, nrow = 256, resolution = 1,
                       center = NULL, decay = 1,
                       invert = FALSE, template = NULL) {
  r <- .chile_raster(ncol, nrow, resolution, template)
  e <- terra::ext(r)
  if (is.null(center)) center <- c((e[1] + e[2]) / 2, (e[3] + e[4]) / 2)
  if (length(center) != 2 || decay <= 0) stop("Invalid center or decay")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  d <- sqrt((xy[, 1] - center[1])^2 + (xy[, 2] - center[2])^2)
  d <- (d / max(d))^decay
  if (invert) d <- 1 - d
  terra::values(r) <- d
  names(r) <- "radial"
  r
}

# -------------------- 10. distance to spatial features -------------------
# Features supplied as two-column XY coordinates, or sampled internally.
nlm_distance <- function(ncol = 256, nrow = 256, resolution = 1,
                         n_features = 10, features = NULL,
                         invert = FALSE, template = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  r <- .chile_raster(ncol, nrow, resolution, template)
  if (is.null(features)) {
    if (n_features < 1) stop("n_features must be positive")
    e <- terra::ext(r)
    features <- cbind(stats::runif(n_features, e[1], e[2]),
                      stats::runif(n_features, e[3], e[4]))
  }
  features <- as.matrix(features)
  if (ncol(features) != 2 || !nrow(features))
    stop("features must be a two-column XY matrix")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  d <- rep(Inf, nrow(xy))
  for (i in seq_len(nrow(features))) {
    q <- (xy[, 1] - features[i, 1])^2 +
         (xy[, 2] - features[i, 2])^2
    d <- pmin(d, q)
  }
  d <- sqrt(d)
  if (invert) d <- max(d) - d
  terra::values(r) <- d
  names(r) <- "distance"
  r
}

# ------------------- 11. marginal distribution transform -----------------
# Rank / Gaussian-quantile mapping preserves rank order, not Pearson
# covariance. Zero inflation introduces ties at zero.
nlm_transform <- function(nlm, distribution = c("normal", "uniform", "lognormal",
                                                  "gamma", "beta", "exponential"),
                          params = NULL, zero_prob = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  distribution <- match.arg(distribution)
  if (zero_prob < 0 || zero_prob > 1) stop("zero_prob must be 0..1")
  r <- nlm[[1]]
  v <- terra::values(r)[, 1]
  ok <- is.finite(v)
  p <- (rank(v[ok], ties.method = "average") - 0.5) / sum(ok)
  # Parameters: use a named list matching the relevant q* function.
  if (is.null(params)) params <- switch(distribution,
    normal = list(mean = 0, sd = 1), uniform = list(min = 0, max = 1),
    lognormal = list(meanlog = 0, sdlog = 1),
    gamma = list(shape = 2, scale = 1),
    beta = list(shape1 = 2, shape2 = 2),
    exponential = list(rate = 1))
  qfun <- switch(distribution,
    normal = stats::qnorm, uniform = stats::qunif,
    lognormal = stats::qlnorm, gamma = stats::qgamma,
    beta = stats::qbeta, exponential = stats::qexp)
  out <- rep(NA_real_, length(v))
  if (length(p)) out[ok] <- do.call(qfun, c(list(p = p), params))
  if (zero_prob > 0) {
    ix <- which(ok)
    out[sample(ix, size = floor(length(ix) * zero_prob))] <- 0
  }
  terra::values(r) <- out
  names(r) <- paste0("transform_", distribution)
  r
}

# ----------------------- 12. reusable field combination -----------------
# Mean and sd standardized within each supplied field; optional weights.
# Allows simple controlled synthetic environmental mixtures.
nlm_mix <- function(..., weights = NULL, mean = 0, sd = 1) {
  fields <- list(...)
  if (!length(fields)) stop("Supply at least one field")
  r <- fields[[1]][[1]]
  for (f in fields[-1]) {
    if (!terra::compareGeom(r, f, stopOnError = FALSE))
      stop("All fields must share a raster grid")
  }
  if (is.null(weights)) weights <- rep(1, length(fields))
  if (length(weights) != length(fields) || all(weights == 0) || sd < 0)
    stop("Invalid weights or sd")
  z <- numeric(terra::ncell(r))
  for (i in seq_along(fields)) {
    x <- terra::values(fields[[i]])[, 1]
    if (any(!is.finite(x))) stop("nlm_mix currently requires complete fields")
    z <- z + weights[i] * as.vector(scale(x))
  }
  z <- mean + sd * as.vector(scale(z))
  terra::values(r) <- z
  names(r) <- "mixture"
  r
}
