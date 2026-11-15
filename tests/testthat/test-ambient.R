test_that("native ambient noise families return valid raster fields", {
  for (kind in c("perlin","simplex","cubic","value","white",
                 "worley","waves","spheres","checkerboard")) {
    x <- nlm_ambient(type=kind,ncol=21,nrow=17,seed=3)
    expect_s4_class(x,"SpatRaster")
    expect_equal(terra::ncell(x),357)
    expect_true(all(is.finite(terra::values(x)[,1])))
  }
})

test_that("ambient anisotropy, warping and worley methods work", {
  x <- nlm_ambient(type="simplex",fractal="fbm",octaves=3,
                   angle=35,anisotropy=4,perturbation="fractal",
                   perturbation_amplitude=9,ncol=30,nrow=25,seed=5)
  y <- nlm_ambient(type="worley",value="distance2sub",distance_ind=c(1,2),
                   ncol=30,nrow=25,seed=5)
  expect_true(all(is.finite(terra::values(x))))
  expect_true(all(is.finite(terra::values(y))))
  expect_equal(terra::values(x),terra::values(nlm_ambient(
    type="simplex",fractal="fbm",octaves=3,angle=35,anisotropy=4,
    perturbation="fractal",perturbation_amplitude=9,ncol=30,nrow=25,seed=5)))
})

test_that("ambient vector fields and galleries are usable", {
  v <- nlm_ambient_flow(ncol=21,nrow=16,seed=1)
  expect_equal(terra::nlyr(v),2)
  expect_identical(names(v),c("dx","dy"))
  grDevices::pdf(file=tempfile(fileext=".pdf"))
  on.exit(grDevices::dev.off(),add=TRUE)
  z <- nlm_ambient_show(type=c("perlin","simplex"),
                        frequency=c(0.01,0.05),ncol=20,nrow=20)
  expect_length(z,4)
})
