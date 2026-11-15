\name{nlm_compare}
\alias{nlm_compare}
\title{Compare observed and simulated scale-dependent statistics}
\description{Computes an observed variogram or K/L/g curve and pointwise simulation quantiles. These are exploratory envelopes, not global Monte Carlo tests.}
\usage{nlm_compare(observed, simulations,
  statistic = c("variogram", "L", "K", "g"), r = NULL,
  lag = NULL, max_dist = NULL, angle = NULL, tolerance = 22.5,
  pairs = 50000, window = NULL, probs = c(0.025, 0.975),
  seed = 1, plot = TRUE)}
\arguments{\item{observed}{Observed single-layer raster or point pattern.}\item{simulations}{Nonempty list of comparable patterns.}\item{statistic}{Variogram, centered L, K, or pair correlation g.}\item{r}{Point-pattern distance scales.}\item{lag, max_dist, angle, tolerance, pairs}{Variogram settings.}\item{window}{Point-pattern rectangular study window.}\item{probs}{Pointwise lower/upper quantiles.}\item{seed}{Seed for sampled variogram pairs.}\item{plot}{Produce a base graphics comparison.}}
\value{Invisibly returns list with curve, observed, simulations and statistic.}
\examples{obs <- nlm_points(ncol = 60, nrow = 60, method = "thomas", n = 60, seed = 1)
sims <- lapply(2:21, function(s) nlm_points(ncol = 60, nrow = 60, method = "thomas", n = 60, seed = s))
cmp <- nlm_compare(obs, sims, statistic = "L")}
