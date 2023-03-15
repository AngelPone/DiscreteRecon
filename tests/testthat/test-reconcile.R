source('./utils.R')



test_that("getVIndexs: get mutable indexs of each row in the A matrix given distance", {
  dts <- prepare_test_data3()
  q <- NROW(dts$hier$incoherent_domain)
  r <- NROW(dts$hier$coherent_domain)
  
  costs <- cal_costeMatrix(dts$hier)
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

optimized_A <- function(hier){
  q <- NROW(hier$incoherent_domain)
  r <- NROW(hier$coherent_domain)
  costs <- cal_costeMatrix(hier)
  indicators <- getVIndexs(costs)
  trained_A <- matrix(0, r, q)
  for (i in seq_along(indicators)){
    trained_A[i, indicators[[i]]] <- rnorm(length(indicators[[i]]))
  }
  trained_A[costs==0] <- 1
  trained_A <- apply(trained_A, 2, function(x){abs(x)/sum(abs(x))})
  trained_A
}

prepare_A <- function(hier){
  r <- hier$quantities['r']
  q <- hier$quantities['q']
  apply(matrix(rnorm(r*q), r, q), 2, function(x){abs(x)/sum(abs(x))})
}

test_that("DFR", {
  dts <- prepare_test_data1()
  basef <- prepare_basef(dts)
  
  A <- dfr(dts$hier, "dfr", 
                 bf_train=basef, obs_train=dts$bts,
                 lambda=0, optimized = FALSE)
  expect_is(A, "dfr")
  expect_equal(colSums(A$A), rep(1, NCOL(A$A)))
  expect_equal(A$method, "dfr")
  
  dts <- prepare_test_data1()
  basef <- prepare_basef(dts)
  foo <- dfr(dts$hier, "dfr", 
                   bf_train=basef, obs_train=dts$bts,
                   optimized = TRUE)
  foo2 <- dfr(dts$hier, "dfr", 
                    bf_train=basef, obs_train=dts$bts,
                    lambda=1, optimized = FALSE)
  # expect_equal(foo$A, foo2$A)
})


test_that("stepwise reconcile: coherentadj",{
  dts <- prepare_test_data1()
  m1 <- structure(list(A = optimized_A(dts$hier), hier=dts$hier,
                       method = "sdfr"), class = "dfr")
  m2 <- structure(list(A = optimized_A(dts$hier), hier=dts$hier,
                       method = "sdfr"), class = "dfr")
  
  prob1 <- reconcile(m1, prepare_basef(dts))
  prob2 <- reconcile(m2, prepare_basef(dts))
  margin_prob2 <- Joint2Marginal(prob2, dts$hier, which = 3)
  adj1 <- coherentadj(prob1, margin_prob2, dts$hier, 3)

  expect_equal(Joint2Marginal(adj1, dts$hier, 3), margin_prob2)
})


test_that("stepwise reconciliation: update_jointDist", {
  s_mat <- rbind(rep(1,2), diag(2))
  lmeta <- dhier(s_mat, matrix(c(0,1,0,2), 2, 2))
  rmeta <- dhier(s_mat, matrix(c(0,1,0,1), 2, 2))
  ldist <- matrix(c(0.1, 0.2, 0.3, 0.2, 0.05, 0.15), 1)
  rdist <- matrix(c(0.3, 0.1, 0.4, 0.2), 1)
  jdist <- update_jointDist(ldist, rdist, lmeta, rmeta)
  expect_equal(jdist,
               matrix(c(0.1, 0.2, 0.06, 0.04, 0.24, 0.16, 0.05, 0.15),1))
})

test_that("stepwise reconcile", {
  dts <- prepare_test_data2(3, 1000)
  
  model <- dfr(dts$hier, method = "sdfr",
               obs_train = dts$bts, bf_train = prepare_basef(dts),
               optimized = TRUE)
  
  expect_is(model, "dfr")
  expect_equal(model$method, "sdfr")
  
  y <- reconcile(model, prepare_basef(dts))
  expect_true(all.equal(rowSums(y), rep(1, 1000)))
})

test_that("topdown reconciliation", {
  
  topdown.train <- function(real){
    coherent_domain <- data.frame(unclass(real$hier$coherent_domain))
    q <- dim(coherent_domain)[1]
    a <- dim(coherent_domain)[2]
    y <- real$bts
    
    concurrence <- vector("numeric", q)
    for (i in 1:dim(y)[1]){
      for (j in 1:q){
        if (all(y[i,] == coherent_domain[j, 2:a])){
          concurrence[j] = concurrence[j] + 1
          break
        }
      }
    }
    tmp <- data.frame(total=coherent_domain$s1, concurrence) %>%
      group_by(total) %>%
      mutate(freq = sum(concurrence), prob = concurrence / freq, n = length(concurrence)) %>%
      mutate(prob = ifelse(is.na(prob), 1/n, prob)) %>%
      pull(prob)
    ds <- sort(unique(coherent_domain[,'s1']))
    A <- matrix(0, q, length(ds))
    colnames(A) <- ds
    for (i in 1:q){
      A[i, as.character(coherent_domain[i, 's1'])] = tmp[i]
    }
    A
  }
  
  
  dts <- prepare_test_data1()
  basef <- prepare_basef(dts)
  r <- dts$hier$quantities['r']
  q <- dts$hier$quantities['q']
  A1 <- topdown.train(dts)
  recf1 <- t(A1 %*% t(basef[[1]]))
  
  A2 <- dfr(dts$hier, method="td", obs_train = dts$bts)
  recf2 <- reconcile(A2, basef)
  expect_equal(Joint2Marginal(recf1, dts$hier, 1), basef[[1]])
  
  expect_equal(recf1, recf2)
})

test_that("BottomUp reconciliation", {
  
  dts <- prepare_test_data1()
  basef <- prepare_basef(dts)
  
  
  A <- dfr(dts$hier, method="bu")
  recf <- reconcile(A, basef)
  expect_equal(Joint2Marginal(recf, dts$hier, 1), 
               marginal2Sum(basef, dts$hier, which = 1))
  expect_equal(Joint2Marginal(recf, dts$hier, 2), 
               basef[[2]])
  expect_equal(Joint2Marginal(recf, dts$hier, 3), 
               basef[[3]])
})

