test_that('directional variogram returns finite scale grid', {
  x <- nlm_anisotropic(ncol=40,nrow=40,seed=1)
  v <- nlm_variogram(x,lag=2,max_dist=12,angle=c(0,90),pairs=4000,seed=1)
  expect_equal(sort(unique(v$angle)),c('0','90'))
  expect_true(all(v$n_pairs>=0))
  expect_true(any(is.finite(v$gamma)))
})

test_that('point diagnostics use the correct study window', {
  p <- nlm_points(ncol=60,nrow=60,n=75,method='binomial',seed=2)
  k <- nlm_k(p,r=c(2,4,8,12))
  expect_equal(nrow(k),4)
  expect_true(all(is.finite(k$K)))
  expect_true(all(diff(k$K)>=0))
})

test_that('patch metrics have sensible areas', {
  p <- nlm_blob(ncol=40,nrow=40,n=5,seed=3)
  m <- nlm_metrics(p)
  expect_true(nrow(m$patches)>0)
  expect_true(all(m$patches$area>0))
  expect_true(all(m$patches$perimeter>0))
})

test_that('simulation comparison produces pointwise bands', {
  obs <- nlm_points(ncol=40,nrow=40,n=35,seed=1)
  sims <- lapply(2:6,function(i)nlm_points(ncol=40,nrow=40,n=35,seed=i))
  cmp <- nlm_compare(obs,sims,statistic='L',r=c(2,4,6),plot=FALSE)
  expect_equal(nrow(cmp$curve),3)
  expect_true(all(cmp$curve$upper>=cmp$curve$lower))
})
