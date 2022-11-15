source('./utils.R')


# source('tests/testthat/utils.R')

# library(CVXR)
# library(Matrix)
# dts <- prepare_test_data2(2, 1000)
# basef <- prepare_basef(dts)
# basef1 <- marginal2Joint(basef, dts$meta, method = "ind")
# q <- dim(dts$meta$incoherent_domain)[1]
# r <- dim(dts$meta$coherent_domain)[1]
# 
# res1 <- allreconcile.train(basef1, dts, optimized = TRUE)
# res <- allreconcile.train(basef1, dts, optimized = FALSE, lambda = 2)
# res2 <- allreconcile.train(basef1, dts, optimized = FALSE, lambda = 1000)
# foo <- reconcile(res, basef, dts$meta)
# foo2 <- reconcile(res2, basef, dts$meta)
# foo1 <- reconcile(res1, basef, dts$meta)
# 
# brier_score(foo1, dts)
# foo2 <- reconcile(res, basef, dts$meta)
# brier_score(foo2, dts)
# costs <- cal_costeMatrix(dts$meta)
# 
# View(unclass(res2))
# indicators <- getVIndexs(costs)
# real_dummy <- cons_realDummy(dts)
# time_window <- 2000
# lambda <- 0
# n <- sum(sapply(indicators, length))
# vs <- Variable(n)
# Dmat <- 2 * construct_Q(r, basef, indicators)
# vec_c <- do.call(c, lapply(1:r, function(i){costs[i,indicators[[i]]]}))
# dvec <- 2 * construct_D(basef, indicators, real_dummy) - time_window * lambda * vec_c
# E_lst <- construct_E(costs, indicators)
# # Dmat <- as.matrix(Dmat)
# 
# Amat <- E_lst$E
# bvec <- E_lst$b0
# A1 <- Amat[1:q, ]
# A2 <- Amat[(q+1):NROW(Amat), ]
# b2 <- bvec[(q+1):NROW(Amat)]
# 
# objective <- Minimize(- t(dvec) %*% vs + 1/2 * quad_form(vs, Dmat))
# 
# constraints <- list(A1 %*% vs == 1, A2 %*% vs >= b2)
# 
# solution <- CVXR::solve(Problem(objective, constraints))
# 
# res2 <- solution$getValue(vs)
# 
# max(abs(t(res1)[t(res1)!=0] - res2[,1]))

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


# test_that("allreconcile.train: matrix form and scalar form",{
#   library(Matrix)
#   dts <- prepare_test_data3()
#   basef <- prepare_basef(dts)
#   
#   cdomain <- dts$meta$coherent_domain
#   idomain <- dts$meta$incoherent_domain
#   
#   costs <- cal_costeMatrix(dts$meta)
#   y <- cons_realDummy(dts)
#   basef <- marginal2Joint(basef, dts$meta, method = "ind")
#   lambda <- 0.001
#   
#   
#   r <- dim(cdomain)[1]
#   q <- dim(idomain)[1]
#   time_window <- dim(y)[1]
#   
#   A <- prepare_A(r, q)
#   approach1 <- function(){
#     Q <- t(basef) %*% basef
#     term1 <- sum(apply(A, 1, function(x){t(x) %*% Q %*% x}))
#     
#     # penalty
#     term2 <- lambda * time_window * sum(sapply(1:r, function(i){
#       sum(costs[i,] * A[i,])
#     }))
#     
#     term3 <- 2 * sapply(1:r, function(i){
#       sum(colSums(basef[y[,i]==1,]) * A[i,])
#     }) %>% sum()
#     
#     term4 <- 1
#     
#     (term1 + term2 - term3 + time_window)/time_window
#   }
#   
#   approach2 <- function(){
#     indicators <- replicate(r, 1:q, simplify = FALSE)
#     vec_c <- as.vector(t(costs))
#     vec_A <- as.vector(t(A))
#     Q <- construct_Q(r, basef)
#     D <- construct_D(basef, indicators, y)
#     vec_z <- as.vector(t(y))
#     
#     term1 <- (t(vec_A) %*% Q %*% vec_A)
#     term2 <- lambda * t(vec_c) %*% vec_A * time_window
#     term3 <- 2 * D %*% vec_A
#     term4 <- t(vec_z) %*% vec_z
#     
#     (term1 + term2 - term3 + term4)/time_window
#   }
#   
#   approach3 <- function(){
#     tmp <- sum((basef %*% t(A) - y)^2) +  lambda * time_window * sum(costs * A)
#     tmp/time_window
#   }
#   
#   expect_equal(approach1(), approach2()[1,1], approach3())
# })

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
  
  # r <- NROW(dts$meta$coherent_domain)
  # q <- NROW(dts$meta$incoherent_domain)
  # costs <- cal_costeMatrix(dts$meta)
  # indicators <- getVIndexs(costs)
  # vec_c <- as.vector(t(costs))
  # A <- optimized_A(dts$meta)
  # time_window <- NROW(basef)
  
  # lambda <- 0.01
  # approach1 <- function(){
  #   Q <- t(basef) %*% basef
  #   term1 <- sum(apply(A, 1, function(x){t(x) %*% Q %*% x}))
  #   
  #   # penalty
  #   term2 <- lambda * time_window * sum(sapply(1:r, function(i){
  #     sum(costs[i,] * A[i,])
  #   }))
  #   
  #   term3 <- 2 * sapply(1:r, function(i){
  #     sum(colSums(basef[y[,i]==1,]) * A[i,])
  #   }) %>% sum()
  #   
  #   term4 <- 1
  #   
  #   (term1 + term2 - term3 + time_window)/time_window
  # }
  # 
  # approach2 <- function(){
  #   vec_A <- do.call(c, sapply(1:r, function(i){A[i, indicators[[i]]]}))
  #   vec_c <- do.call(c, sapply(1:r, function(i){costs[i, indicators[[i]]]}))
  #   Q <- construct_Q(r, basef, indicators)
  #   zd <- 2 * construct_D(basef, indicators, y)
  #   lc <- lambda * time_window * vec_c
  #   
  #   term1 <- t(vec_A) %*% Q %*% vec_A
  #   term2 <- sum(lc * vec_A)
  #   term3 <- sum(zd * vec_A)
  #   term4 <- time_window
  #   
  #   ((term1 + term2 - term3 + term4)/time_window)[1,1]
  # }
  # 
  # expect_equal(approach1(), approach2())
  
  # vec_A <- t(A)[t(A)>0]
  # # construct_E
  # E <- construct_E(costs, indicators)$E[1:q,]
  # 
  # expect_equal(as.numeric(E %*% vec_A), rep(1, q))
  # 
  # # construct_res
  # vec_A <- t(A)[t(A)>0]
  # expect_equal(construct_res(costs, vec_A, indicators), A)
  
  # workflow
  foo <- allreconcile.train(basef, dts, optimized=TRUE)
  foo2 <- allreconcile.train(basef, dts, 100, FALSE)
  expect_equal(foo, foo2)
})


test_that("stepwise reconcile: coherentadj",{
  dts <- prepare_test_data1()
  prob1 <- reconcile(optimized_A(dts$meta), prepare_basef(dts), dts$meta)
  prob2 <- reconcile(optimized_A(dts$meta), prepare_basef(dts), dts$meta)
  margin_prob2 <- Joint2Marginal(prob2, dts$meta, which = 3)
  adj1 <- coherentadj(prob1, margin_prob2, dts$meta, 3)
  
  expect_equal(Joint2Marginal(adj1, dts$meta, 3), margin_prob2)
})

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

