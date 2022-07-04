# test dhts object and domain

test_that("domain", {
  bts <- matrix(sample(0:1, 100, replace=TRUE), ncol=2)
  domain <- matrix(c(0, 1, 0, 1), 2)
  s_mat <- matrix(c(1, 1, 0, 1, 0, 1), 3)
  a <- dhts(bts, s_mat, domain)
  t1 <- a$domain$coherent_domain
  t2 <- a$domain$incoherent_domain
  
  expect_is(t1, "coherent_domain")
  expect_is(t2, "incoherent_domain")
  expect_equal(dim(t1), c(4,3))
  expect_equal(dim(t2), c(12, 3))
  expect_equal(t1[,1], t1[,2] + t1[,3])
  expect_equal(colnames(t1), paste0("s", 1:3))
  expect_equal(colnames(t2), paste0("s", 1:3))
  expect_equal(max(t1), 2)
  expect_equal(min(t1), 0)
})



test_that("domain", {
  bts <- matrix(sample(0:1, 700, replace=TRUE), ncol=7)
  domain <- rbind(rep(0, 7), rep(1, 7))
  s_mat <- rbind(rep(1,7), diag(7))
  a <- dhts(bts, s_mat, domain)
  t1 <- a$domain$coherent_domain
  t2 <- a$domain$incoherent_domain
  expect_is(t1, "coherent_domain")
  expect_is(t2, "incoherent_domain")
  expect_equal(dim(t1), c(128, 8))
  expect_equal(dim(t2), c(128*8, 8))
  expect_equal(t1[,1], rowSums(t1[,2:8]))
  expect_equal(colnames(t1), paste0("s", 1:8))
  expect_equal(colnames(t2), paste0("s", 1:8))
  expect_equal(max(t1), 7)
  expect_equal(min(t1), 0)
})


# different domain for bottom level time series

test_that("domain", {
  upper <- sample(3:8, 7, replace = TRUE)
  domain <- rbind(upper, upper + sample(1:2, 7, replace = TRUE))
  s_mat <- rbind(rep(1,7), diag(7))
  bts <- apply(domain, 2, function(x){
    sample(x[1]:x[2], size=100, replace = TRUE)
  })
  a <- dhts(bts, s_mat, domain)
  t1 <- a$domain$coherent_domain
  t2 <- a$domain$incoherent_domain
  
  dslower <- apply(domain, 2, function(x){(x[2] - x[1] + 1)})
  dsranges <- (s_mat %*% t(domain))[,2] - (s_mat %*% t(domain))[,1] + 1
  expect_is(t1, "coherent_domain")
  expect_is(t2, "incoherent_domain")
  expect_equal(dim(t1), c(prod(dslower), 8))
  expect_equal(dim(t2), as.numeric(c(prod(dsranges), 8)))
  expect_equal(t1[,1], rowSums(t1[,2:8]))
  expect_equal(colnames(t1), paste0("s", 1:8))
  expect_equal(colnames(t2), paste0("s", 1:8))
  expect_equal(max(t1), as.numeric(rowSums(domain)[2]))
  expect_equal(min(t1), as.numeric(min(domain)))
})



test_that("multiple levels example", {
  
  upper <- sample(3:8, 7, replace = TRUE)
  domain <- rbind(upper, upper + sample(1:2, 7, replace = TRUE))
  s_mat <- rbind(c(1,1,1,1,1,1,1), c(1,1,1,0,0,0,0), c(0,0,0,1,1,1,1), diag(7))
  bts <- apply(domain, 2, function(x){
    sample(x[1]:x[2], size=100, replace = TRUE)
  })
  dslower <- apply(domain, 2, function(x){(x[2] - x[1] + 1)})
  dsranges <- (s_mat %*% t(domain))[,2] - (s_mat %*% t(domain))[,1] + 1
  a <- dhts(bts, s_mat, domain)
  t1 <- a$domain$coherent_domain
  t2 <- a$domain$incoherent_domain
  
  expect_equal(dim(t1), c(prod(dslower), 10))
  expect_equal(dim(t2), as.numeric(c(prod(dsranges), 10)))
  expect_equal(t1[,1], rowSums(t1[,4:10]))
  expect_equal(t1[,2], rowSums(t1[,4:6]))
  expect_equal(t1[,3], rowSums(t1[,7:10]))
  expect_equal(colnames(t1), paste0("s", 1:10))
  expect_equal(colnames(t2), paste0("s", 1:10))
  expect_equal(max(t1), as.numeric(rowSums(domain)[2]))
  expect_equal(min(t1), as.numeric(min(domain)))
  expect_equal(attr(t1, 'm'), 7)
  expect_equal(attr(t2, 'm'), 7)
})

