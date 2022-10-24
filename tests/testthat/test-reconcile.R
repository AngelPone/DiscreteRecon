source('./utils.R')

test_that("getVIndexs: get mutable indexs of each row in the A matrix given distance", {
  dts <- prepare_test_data3()
  q <- dim(dts$meta$incoherent_domain)[1]
  r <- dim(dts$meta$coherent_domain)[1]
  
  costs <- cal_costeMatrix(dts$meta)
  vIndexs <- getVIndexs(costs)
  expect_equal(length(vIndexs), r)
  
  mins <- list()
  for (i in 1:q){
    mins[[i]] <- which((costs[,i] == min(costs[,i])))
  }
  for (i in seq_along(vIndexs)){
    included <- (1:q)[sapply(mins, function(x){i %in% x})]
    expect_equal(length(setdiff(vIndexs[[i]], included)),0)
  }
})

optimized_A <- function(meta){
  q <- NROW(dts$meta$incoherent_domain)
  r <- NROW(dts$meta$coherent_domain)
  costs <- cal_costeMatrix(meta)
  indicators <- getVIndexs(costs)
  trained_A <- matrix(0, r, q)
  for (i in seq_along(indicators)){
    trained_A[i, indicators[[i]]] <- rnorm(length(indicators[[i]]))
  }
  trained_A[costs==0] <- 1
  trained_A <- apply(trained_A, 2, function(x){abs(x)/sum(abs(x))})
  structure(trained_A, class="rec_mat")
}

prepare_A <- function(r, q){
  apply(matrix(rnorm(r*q), r, q), 2, function(x){abs(x)/sum(abs(x))})
}


test_that("allreconcile.train: matrix form and scalar form",{
  library(Matrix)
  dts <- prepare_test_data3()
  basef <- prepare_basef(dts)
  
  cdomain <- dts$meta$coherent_domain
  idomain <- dts$meta$incoherent_domain
  
  costs <- cal_costeMatrix(dts$meta)
  y <- cons_realDummy(dts)
  basef <- marginal2Joint(basef, dts$meta, method = "ind")
  lambda <- 0.001
  
  
  r <- dim(cdomain)[1]
  q <- dim(idomain)[1]
  time_window <- dim(y)[1]
  
  A <- prepare_A(r, q)
  approach1 <- function(){
    Q <- t(basef) %*% basef
    term1 <- sum(apply(A, 1, function(x){t(x) %*% Q %*% x}))
    
    # penalty
    term2 <- lambda * time_window * sum(sapply(1:r, function(i){
      sum(costs[i,] * A[i,])
    }))
    
    term3 <- 2 * sapply(1:r, function(i){
      sum(colSums(basef[y[,i]==1,]) * A[i,])
    }) %>% sum()
    
    term4 <- 1
    
    (term1 + term2 - term3 + time_window)/time_window
  }
  
  approach2 <- function(){
    indicators <- replicate(r, 1:q, simplify = FALSE)
    vec_c <- as.vector(t(costs))
    vec_A <- as.vector(t(A))
    Q <- construct_Q(r, basef)
    D <- construct_D(basef, indicators, y)
    vec_z <- as.vector(t(y))
    
    term1 <- (t(vec_A) %*% Q %*% vec_A)
    term2 <- lambda * t(vec_c) %*% vec_A * time_window
    term3 <- 2 * D %*% vec_A
    term4 <- t(vec_z) %*% vec_z
    
    (term1 + term2 - term3 + term4)/time_window
  }
  
  approach3 <- function(){
    tmp <- sum((basef %*% t(A) - y)^2) +  lambda * time_window * sum(costs * A)
    tmp/time_window
  }
  
  
  # indicators <- getVIndexs(costs)
  
  expect_equal(approach1(), approach2()[1,1], approach3())
  
  
  # Dmat <- construct_Q(r, basef, indicators)
  # D <- construct_D(basef, indicators, y)
  # dvec <- t(as.matrix(2 * D))
  # vec_z <- matrix(t(y), ncol = 1)
  # trained_A <- optimized_A(r, q, indicators, costs)
  # vec_A <- matrix(t(trained_A), ncol = 1)
  # 
  # #检验转换正确
  # expect_equal(vec_A[1:720,1], trained_A[1,])
  # 
  # # 检验矩阵形式的正确性
  # mat_form <- 1/time_window * t(vec_A) %*% construct_Q(r, basef) %*% vec_A +
  #   (lambda * vec_c - 2/time_window * D) %*% vec_A +
  #   1/time_window * t(vec_z) %*% vec_z
  # mat_form <- mat_form[1,1]
  # 
  # scalar_form <- sum((t(trained_A %*% t(basef)) - y)^2)/time_window + sum(costs * trained_A)
  # 
  # expect_equal(scalar_form, mat_form)
})

test_that("allreconcile.train: optimize=FALSE", {
  dts <- prepare_test_data1()
  basef <- prepare_basef(dts)
  basef <- marginal2Joint(basef, dts$meta, method = "ind")
  
  A <- allreconcile.train(basef, dts, lambda=0, optimized = FALSE)
  expect_is(A, "rec_mat")
  expect_equal(colSums(A), rep(1, NCOL(A)))
})

test_that("allreconcile.train: optimize=TRUE", {
  dts <- prepare_test_data2(4)
  y <- cons_realDummy(dts)
  basef <- prepare_basef(dts)
  basef <- marginal2Joint(basef, dts$meta, method = "ind")
  
  r <- NROW(dts$meta$coherent_domain)
  q <- NROW(dts$meta$incoherent_domain)
  costs <- cal_costeMatrix(dts$meta)
  indicators <- getVIndexs(costs)
  vec_c <- as.vector(t(costs))
  A <- optimized_A(r, q, indicators, costs)
  time_window <- NROW(basef)
  
  lambda <- 0.01
  approach1 <- function(){
    Q <- t(basef) %*% basef
    term1 <- sum(apply(A, 1, function(x){t(x) %*% Q %*% x}))
    
    # penalty
    term2 <- lambda * time_window * sum(sapply(1:r, function(i){
      sum(costs[i,] * A[i,])
    }))
    
    term3 <- 2 * sapply(1:r, function(i){
      sum(colSums(basef[y[,i]==1,]) * A[i,])
    }) %>% sum()
    
    term4 <- 1
    
    (term1 + term2 - term3 + time_window)/time_window
  }
  
  approach2 <- function(){
    vec_A <- do.call(c, sapply(1:r, function(i){A[i, indicators[[i]]]}))
    vec_c <- do.call(c, sapply(1:r, function(i){costs[i, indicators[[i]]]}))
    Q <- construct_Q(r, basef, indicators)
    zd <- 2 * construct_D(basef, indicators, y)
    lc <- lambda * time_window * vec_c
    
    term1 <- t(vec_A) %*% Q %*% vec_A
    term2 <- sum(lc * vec_A)
    term3 <- sum(zd * vec_A)
    term4 <- time_window
    
    ((term1 + term2 - term3 + term4)/time_window)[1,1]
  }
  
  expect_equal(approach1(), approach2())
  
  vec_A <- t(A)[t(A)>0]
  # construct_E
  E <- construct_E(costs, indicators)$E[1:q,]
  
  expect_equal(as.numeric(E %*% vec_A), rep(1, q))
  
  # construct_res
  vec_A <- t(A)[t(A)>0]
  expect_equal(construct_res(costs, vec_A, indicators), A)
  
  # workflow
  foo <- allreconcile.train(basef, dts, optimized=TRUE)
  foo2 <- allreconcile.train(basef, dts, 100, FALSE)
  expect_equal(foo, foo2)
})


test_that("stepwise reconcile: coherentadj"){
  dts <- prepare_test_data1()
  prob1 <- reconcile(optimized_A(dts$meta), prepare_basef(dts), dts$meta)
  prob2 <- reconcile(optimized_A(dts$meta), prepare_basef(dts), dts$meta)
  margin_prob2 <- Joint2Marginal(prob2, dts$meta, which = 3)
  adj1 <- coherentadj(prob1, margin_prob2, dts$meta, 3)
  
  expect_equal(Joint2Marginal(adj1, dts$meta, 3), margin_prob2)
}

test_that("stepwise reconcile", {
  dts <- prepare_test_data2(3, 1000)
  basef <- prepare_basef(dts)
  
  model <- reconcile_train(basef, dts, step_wise = TRUE, optimized = FALSE, lambda = 1)
  
  expect_is(model, "sw_res")
  
  y <- reconcile(model, basef)
  expect_true("coherent" %in% class(y))
  expect_true("jdist" %in% class(y))
  expect_true("rec" %in% class(y))
  expect_equal(NROW(dts$meta$coherent_domain), NCOL(y))
})

test_that("topdown reconciliation", {
  dts <- prepare_test_data2(4)
  r = dim(dts$meta$coherent_domain)[1]
  q <- length(unique(dts$meta$incoherent_domain[,'s1']))
  A <- topdown.train(dts)
  expect_equal(dim(A), c(r, q))
  expect(all(!is.na(A)), "shoule be na")
  
  basef <- prepare_basef(dts)
  recf <- reconcile(A, basef)
  expect_true("coherent" %in% class(recf))
  expect_true("jdist" %in% class(recf))
  expect_true("td" %in% class(recf))
  expect_equal(Joint2Marginal(recf, dts$meta, 1), basef[[1]])
})

