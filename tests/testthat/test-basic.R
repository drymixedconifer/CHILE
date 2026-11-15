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

test_that('revised field families generate finite realizations', {
  for (f in list(
    nlm_gaussianfield(ncol=33,nrow=25,range=5,seed=1),
    nlm_midpoint(ncol=33,nrow=25,seed=1),
    nlm_ambient(ncol=33,nrow=25,type='perlin',seed=1),
    nlm_ambient(ncol=33,nrow=25,type='simplex',seed=1),
    nlm_ambient(ncol=33,nrow=25,type='value',seed=1),
    nlm_ambient(ncol=33,nrow=25,type='worley',seed=1),
    nlm_ambient(ncol=33,nrow=25,type='white',seed=1)
  )) {
    expect_s4_class(f,'SpatRaster')
    expect_equal(terra::ncell(f),825)
    expect_true(all(is.finite(terra::values(f)[,1])))
  }
})

test_that('field marginal parameters and distinct noise types behave', {
  x <- nlm_gaussianfield(ncol=40,nrow=30,range=5,mean=12,sill=9,seed=2)
  expect_equal(mean(terra::values(x)[,1]),12,tolerance=1e-6)
  expect_equal(stats::sd(terra::values(x)[,1]),3,tolerance=1e-6)
  a <- nlm_ambient(ncol=40,nrow=30,type='perlin',seed=3)
  b <- nlm_ambient(ncol=40,nrow=30,type='value',seed=3)
  c <- nlm_ambient(ncol=40,nrow=30,type='simplex',seed=3)
  expect_false(isTRUE(all.equal(terra::values(a),terra::values(b))))
  expect_false(isTRUE(all.equal(terra::values(b),terra::values(c))))
})

test_that('default CHILE gallery renders without error', {
  grDevices::pdf(file=tempfile(fileext='.pdf'))
  on.exit(grDevices::dev.off(),add=TRUE)
  expect_silent(nlm_show(ncol=32,nrow=32))
})
