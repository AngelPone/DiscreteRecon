source('./utils.R')

test_that("getVIndexs: get mutable indexs of each row in the A matrix given distance", {
  dts <- prepare_test_data3()
  q <- dim(dts$domain$incoherent_domain)[1]
  r <- dim(dts$domain$coherent_domain)[1]
  
  costs <- cal_costeMatrix(dts$domain$incoherent_domain, dts$domain$coherent_domain)
  vIndexs <- getVIndexs(costs)
  expect_equal(length(vIndexs), r)
  
  mins <- list()
  for (i in 1:q){
    mins[[i]] <- which((costs[,i] == min(costs[,i])) & (min(costs[,i]) > 0))
  }
  for (i in seq_along(vIndexs)){
    included <- (1:q)[sapply(mins, function(x){i %in% x})]
    expect_equal(length(setdiff(vIndexs[[i]], included)),0)
  }
})

optimized_A <- function(r, q, indicators, costs){
  trained_A <- matrix(0, r, q)
  for (i in seq_along(indicators)){
    trained_A[i, indicators[[i]]] <- rnorm(length(indicators[[i]]))
  }
  trained_A[costs==0] <- 1
  trained_A <- apply(trained_A, 2, function(x){abs(x)/sum(abs(x))})
  trained_A
}


test_that("allreconcile.train: matrix form and scalar form",{
  dts <- prepare_test_data3()
  basef <- prepare_basef(dts)
  
  cdomain <- dts$domain$coherent_domain
  idomain <- dts$domain$incoherent_domain
  
  costs <- cal_costeMatrix(idomain, cdomain)
  y <- cons_realDummy(dts)
  basef <- marginal2Joint(basef, idomain)
  
  r <- dim(cdomain)[1]
  q <- dim(idomain)[1]
  time_window <- dim(y)[1]
  lambda <- 1
  
  vec_c <- matrix(t(costs), nrow = 1)
  
  # indicators <- getVIndexs(costs)
  indicators <- replicate(r, 1:q, simplify = FALSE)
  
  Dmat <- 2 * construct_Q(r, basef, indicators)
  D <- construct_D(basef, indicators, y)
  dvec <- t(as.matrix(2 * D))
  vec_z <- matrix(t(y), ncol = 1)
  trained_A <- optimized_A(r, q, indicators, costs)
  vec_A <- matrix(t(trained_A), ncol = 1)
  
  #检验转换正确
  expect_equal(vec_A[1:720,1], trained_A[1,])
  
  # 检验矩阵形式的正确性
  mat_form <- 1/time_window * t(vec_A) %*% construct_Q(r, basef) %*% vec_A +
    (lambda * vec_c - 2/time_window * D) %*% vec_A +
    1/time_window * t(vec_z) %*% vec_z
  mat_form <- mat_form[1,1]
  
  scalar_form <- sum((t(trained_A %*% t(basef)) - y)^2)/time_window + sum(costs * trained_A)
  
  expect_equal(scalar_form, mat_form)
  
})

test_that("allreconcile.train: optimized matrix form and scalar form",{
  dts <- prepare_test_data3()
  basef <- prepare_basef(dts)
  
  cdomain <- dts$domain$coherent_domain
  idomain <- dts$domain$incoherent_domain
  
  costs <- cal_costeMatrix(idomain, cdomain)
  y <- cons_realDummy(dts)
  basef <- marginal2Joint(basef, idomain)
  
  r <- dim(cdomain)[1]
  q <- dim(idomain)[1]
  time_window <- dim(y)[1]
  lambda <- 1
  
  
  indicators <- getVIndexs(costs)
  
  
  for (i in 1:r){
    indicators[[i]] <- sort(c(indicators[[i]], which(costs[i,] == 0)))
  }
  
  vec_c <- NULL
  for (i in 1:length(indicators)){
    vec_c <- c(vec_c, costs[i,indicators[[i]]])
  }
  vec_c <- matrix(vec_c, ncol = 1)
  
  
  Dmat <- 2 * construct_Q(r, basef, indicators)
  D <- construct_D(basef, indicators, y)
  dvec <- t(as.matrix(2 * D))
  
  trained_A <- optimized_A(r, q, indicators, costs)
  vec_A <-  matrix(t(trained_A)[t(trained_A)>0], ncol = 1)
  
  # 检验矩阵形式的正确性
  mat_form <- 1/time_window * t(vec_A) %*% construct_Q(r, basef, indicators) %*% vec_A +
    (t(lambda * vec_c) - 2/time_window * D) %*% vec_A + 1
  mat_form <- mat_form[1,1]
  
  scalar_form <- sum((t(trained_A %*% t(basef)) - y)^2)/time_window + sum(costs * trained_A)
  
  expect_equal(scalar_form, mat_form)
  
})


test_that("all reconcile", {
  dts <- prepare_test_data1()
  basef <- prepare_basef(dts)
  
  model <- reconcile_train(basef, dts, step_wise = FALSE)
  
  expect_is(model, "rec_mat")
  expect_equal(colSums(model$A), rep(1, 12))
  
  y <- reconcile(model, basef)
  expect_is(y, "jdist-rec")
  expect_equal(dim(y), c(100, 4))
})

test_that("stepwise reconcile", {
  dts <- prepare_test_data2(4)
  basef <- prepare_basef(dts)
  
  model <- reconcile_train(basef, dts)
  
  expect_is(model, "sw_res")
  
  y <- reconcile(model, basef)
  expect_is(y, 'jdist-rec')
})

test_that("topdown reconciliation", {
  dts <- prepare_test_data2(4)
  r = dim(dts$domain$coherent_domain)[1]
  q <- length(unique(dts$domain$incoherent_domain[,'s1']))
  A <- topdown.train(dts)
  expect_equal(dim(A), c(r, q))
  expect(all(!is.na(A)), "shoule be na")
  
  basef <- prepare_basef(dts)
  recf <- reconcile(A, basef)
  expect_is(recf, 'jdist-td')
  expect_equal(Joint2Marginal(recf, dts$domain$coherent_domain, 1), basef[[1]])
})

