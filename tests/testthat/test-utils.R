source('./utils.R')

test_that("convert dhts to dummies", {
  sample_dhts <- prepare_test_data1()
  domain <- sample_dhts$hier$coherent_domain
  m <- dim(sample_dhts$hier$domain)[2]
  n <- dim(domain)[2]
  a <- cons_realDummy(sample_dhts$bts, sample_dhts$hier)
  for (i in 1:NROW(a)){
    expect_equal(as.numeric(domain[which.max(a[i,]), (n-m+1):n]), 
                 as.numeric(sample_dhts$bts[i,]))
  }
})

test_that("calculate cost matrix", {
  sample_dhts <- prepare_test_data1()
  
  d2 <- sample_dhts$hier$coherent_domain
  d1 <- sample_dhts$hier$incoherent_domain
  distance <- cal_costeMatrix(sample_dhts$hier)
  
  expect_equal(dim(distance), c(dim(d2)[1], dim(d1)[1]))
  for (i in 1:dim(d2)[1]){
    for (j in 1:dim(d1)[1]){
      expect_equal(distance[i, j], sum(abs(d2[i,] - d1[j,])))
    }
  }
  expect_equal(apply(distance[,sample_dhts$hier$coherent_flags$coherent], 1, min), rep(0, NROW(d2)))
})

test_that("Marginal distributions to Joint distributions", {
  sample_dhts <- prepare_test_data1()
  basef <- prepare_basef(sample_dhts)
  r <- dim(sample_dhts$hier$coherent_domain)[1]
  q <- dim(sample_dhts$hier$incoherent_domain)[1]
  n <- dim(sample_dhts$hier$coherent_domain)[2]
  m <- dim(sample_dhts$hier$domain)[2]
  
  f1 <- marginal2Joint(basef, sample_dhts$hier, method = "ind")
  expect_equal(dim(f1), c(100, q))
  for (i in 1:q){
    a <- sample_dhts$hier$incoherent_domain[i,] %>% as.numeric()
    b <- NULL
    for (j in 1:length(basef)){
      b <- cbind(b, basef[[j]][, as.character(a[j])])
    }
    expect_equal(
      f1[,i],
      apply(b, 1, prod)
    )
  }
  
  f1 <- marginal2Joint(basef, sample_dhts$hier, method = "bu")
  for (i in 1:r){
    a <- sample_dhts$hier$coherent_domain[i,(n-m+1):n] %>% as.numeric()
    b <- NULL
    for (j in (n-m+1):(n)){
      b <- cbind(b, basef[[j]][, as.character(a[j - (n-m)])])
    }
    expect_equal(
      f1[,i],
      apply(b, 1, prod)
    )
  }
  expect_equal(dim(f1), c(100, r))
})

test_that("Joint Distributions to Marginal distribution of specific dimension for bottom up distribution", {
  dts <- prepare_test_data1()
  basef <- prepare_basef(dts)
  domain <- dts$hier$coherent_domain
  
  n <- NROW(dts$hier$s_mat)
  m <- NCOL(dts$hier$s_mat)
  jdist <- marginal2Joint(basef, dts$hier, method = "bu")
  mdist <- Joint2Marginal(jdist, dts$hier)
  for (i in (n-m+1):length(mdist)){
    expect_equal(mdist[[i]], basef[[i]])
    expect_equal(colnames(mdist[[i]]), colnames(basef[[i]]))
  }
})

test_that("Joint Distributions to Marginal distribution of specific dimension for reconciled distribution", {
  dts <- prepare_test_data1()
  recdist <- prepare_recdist(dts)

  mdist <- Joint2Marginal(recdist, dts$hier)
  
  domain <- dts$hier$coherent_domain
  for(i in seq_along(mdist)){
    cdist <- mdist[[i]]
    
    for (col in colnames(cdist)){
      expect_equal(rowSums(recdist[,which(domain[,i] == col),drop=FALSE]),
                   cdist[,col])
    }
  }
})


test_that("Marginal distributions to sum assuming independence", {
  dts <- prepare_test_data3()
  basef <- prepare_basef(dts)
  domain <- dts$hier$coherent_domain
  res <- marginal2Sum(basef, dts$hier)
  
  expect_equal(length(res), 7)
  expect_equal(dim(res[[2]]), c(100, 3))
  expect_equal(dim(res[[3]]), c(100, 3))
  expect_equal(dim(res[[1]]), c(100, 5))
  
  expect_equal(res[[1]][,"0"], apply(
    do.call(cbind, lapply(basef[4:7], function(x){x[,"0"]})), 1, prod
  ))
  expect_equal(res[[1]][,"4"], apply(
    do.call(cbind, lapply(basef[4:7], function(x){x[,"1"]})), 1, prod
  ))
  
  res <- marginal2Sum(basef, dts$hier, 1)
  expect_equal(dim(res), c(100, 5))
  
  res <- marginal2Sum(basef, dts$hier, 2)
  expect_equal(dim(res), c(100, 3))
  
  res <- marginal2Sum(basef, dts$hier, 3)
  expect_equal(dim(res), c(100, 3))
  expect_equal(res[,"2"], apply(
    do.call(cbind, lapply(basef[6:7], function(x){x[,"1"]})), 1, prod
  ))
  expect_equal(res[,"0"], apply(
    do.call(cbind, lapply(basef[6:7], function(x){x[,"0"]})), 1, prod
  ))
  expect_equal(res[,"1"], 
                basef[[6]][,"0"] * basef[[7]][,"1"] + 
                basef[[6]][,"1"] * basef[[7]][,"0"])
})





