source('./utils.R')

test_that("getCoherentFlag", {
  dts <- prepare_test_data1()
  flags <- getCoherentFlag(dts$domain$incoherent_domain, dts$domain$coherent_domain)
  expect_equal(flags$incoherent, c(2,3,4,6, 7, 9,10,11))
  expect_equal(flags$coherent, c(1, 5, 8, 12))
})

test_that("brier Score", {
  dts <- prepare_test_data1()
  jdist_rec <- prepare_recdist(dts)
  brier_score(jdist_rec, dts)
  
  basef <- marginal2Joint(prepare_basef(dts), dts$domain$incoherent_domain)
  brier_score(basef, dts)
})

test_that("point forecasts", {
  dts <- prepare_test_data1()
  jdist_rec <- prepare_recdist(dts)
  
  pointf <- point_forecast(jdist_rec, dts$domain$coherent_domain)
  expect_equal(pointf[,1], pointf[,2] + pointf[,3])
  
  mae <- function(x, y){mean(abs(x-y))}
  rmse <- function(x, y){mean((x-y)^2)}
  
  point_metric(jdist_rec, dts, mae)
  point_metric(jdist_rec, dts, rmse)
})