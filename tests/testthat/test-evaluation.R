source('./utils.R')

test_that("getCoherentFlag", {
  dts <- prepare_test_data1()
  flags <- getCoherentFlag(dts$hier$incoherent_domain, dts$hier$coherent_domain, dts$hier$s_mat)
  expect_equal(flags$incoherent, c(2,3,4,6, 7, 9,10,11))
  expect_equal(flags$coherent, c(1, 5, 8, 12))
})

test_that("brier Score", {
  dts <- prepare_test_data1()
  
  jdist_rec <- prepare_recdist(dts)
  brier_score(jdist_rec, dts$bts, dts$hier)
  
  # base
  brier_score(prepare_basef(dts), dts$bts, dts$hier)
})

test_that("point forecasts", {
  dts <- prepare_test_data1()
  jdist_rec <- prepare_recdist(dts)
  
  pointf <- point_forecast(jdist_rec, dts$hier)
  expect_equal(pointf[,1], pointf[,2] + pointf[,3])
  
  mae <- function(x, y){mean(abs(x-y))}
  rmse <- function(x, y){mean((x-y)^2)}
  # point metric for reconciled forecasts
  point_metric(jdist_rec, dts$bts, dts$hier, mae)
  point_metric(jdist_rec, dts$bts, dts$hier, rmse)
  
  # bottom up
  basef <- prepare_basef(dts)
  bu_mae <- point_metric(marginal2Joint(basef, dts$hier, method="bu"), dts$bts, 
                         dts$hier, mae)
  
  # basef
  basef_mae <- point_metric(basef, dts$bts, dts$hier,mae)
  
  expect_equal(bu_mae[2:3], basef_mae[2:3])
})
