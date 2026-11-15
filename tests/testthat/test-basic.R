test_that('generators run without external environment', {
  x <- nlm_anisotropic(ncol=32,nrow=24,seed=1)
  expect_s4_class(x,'SpatRaster')
  expect_equal(terra::ncell(x),768)
  expect_equal(terra::nlyr(nlm_crosscorrelated(ncol=24,nrow=24)),2)
  expect_s4_class(nlm_polygon(ncol=32,nrow=24,n=5)$raster,'SpatRaster')
  expect_s4_class(nlm_blob(ncol=32,nrow=24,n=5)$raster,'SpatRaster')
  expect_s4_class(nlm_filament(ncol=32,nrow=24,n=5)$raster,'SpatRaster')
  expect_s4_class(nlm_partition(ncol=32,nrow=24)$raster,'SpatRaster')
})
test_that('point process returns coordinates and count raster', {
  p <- nlm_points(ncol=35,nrow=30,method='binomial',n=60,seed=1)
  expect_equal(nrow(p$points),60)
  expect_equal(sum(terra::values(p$raster)),60)
  for(m in c('thomas','matern_cluster','inhibition','strauss','lgcp',
             'grid','jittered','anisotropic_cluster','mixture')) {
    z <- nlm_points(ncol=40,nrow=40,n=30,method=m,seed=1)
    expect_s4_class(z$raster,'SpatRaster')
    expect_equal(sum(terra::values(z$raster)),nrow(z$points))
  }
})
