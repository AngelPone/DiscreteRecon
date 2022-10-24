#' reconcile base probabilistic forecasts
#' 
#' Method for reconciling base probabilistic forecasts 
#' @import Matrix
#' @param basef incoherent joint base probabilistic forecasts
#' @param real dhts object holding observed series
#' @param lambda penalty terms of cost, default is 1.
#' @param optimized indicating whether optimized. i.e., only moving probabilities to nearest coherent points.
#' @tag multiple-levels
#' @return r*q transformation matrix with class `rec_mat` 
allreconcile.train <- function(basef, 
                               real,
                               lambda = 0,
                               optimized=TRUE) {
  stopifnot(is(real, "dhts"))
  stopifnot("incoherent" %in% class(basef))
  stopifnot("jdist" %in% class(basef))
  
  real_dummy <- cons_realDummy(real)
  stopifnot(NROW(basef) == NROW(real_dummy))
  
  cdomain <- real$meta$coherent_domain
  idomain <- real$meta$incoherent_domain
  
  r <- NROW(cdomain)
  q <- NROW(idomain)
  time_window <- NROW(real_dummy)
  costs <- cal_costeMatrix(real$meta)
  
  vec_c <- as.vector(t(costs))
  
  if (optimized){
    lambda <- 0
    indicators <- getVIndexs(costs)
    Dmat <- 2 * construct_Q(r, basef, indicators)
    vec_c <- do.call(c, lapply(1:r, function(i){costs[i,indicators[[i]]]}))
    dvec <- 2 * construct_D(basef, indicators, real_dummy) - 
      time_window * lambda * vec_c
    
    E_lst <- construct_E(costs, indicators)
    
    Amat <- t(E_lst$E)
    bvec <- E_lst$b0
    n_eq <- q
  } else {
    Dmat <- 2 * construct_Q(r, basef)  
    dvec <- 2 * construct_D(basef, NULL, real_dummy) - time_window * lambda * vec_c
    
    A1 <- matrix(replicate(r, diag(q)), q)
    b1 <- rep(1, q)
    
    A2 <- diag(r * q)
    b2 <- rep(0, r * q)
    
    A3 <- -diag(r * q)
    b3 <- rep(-1, r * q)
    Amat <- t(rbind(A1, A2, A3))
    bvec <- c(b1, b2, b3)
    n_eq <- q
  }

  solution <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = n_eq)$solution
  solution[solution<1e-10] = 0
  if (optimized){
    Ahat <- construct_res(costs, solution, indicators)
  }else{
    Ahat <- matrix(solution, r, q, byrow = TRUE)
  }
  structure(apply(Ahat, 2, function(x){x/sum(x)}), class="rec_mat")
}


#' function to train reconciliation model
#' 
#' @param basef list containing base probabilistic forecasts of each series
#' @param real dhts object containing corresponding observations of basef
#' @param lambda penality when optmized=FALSE
#' @param optimized boolean indicating if optimization 
#' @param step_wise boolean indicating if reconciliation step wisely.
#' @tag single-level
#' @return sw_res object or rec_mat object
#' @export
reconcile_train <- function(basef, real, 
                            lambda=0, 
                            optimized=TRUE, 
                            step_wise=TRUE){
  if (optimized) lambda <- 0 
  if (!step_wise) {
    basef <- marginal2Joint(basef, real$meta, method="ind")
    return(allreconcile.train(basef, real, lambda=lambda, 
                              optimized = optimized))
  } else {
    return(stepwise_reconcile.train(basef, real, lambda = lambda, optimized = optimized))
  }
}

#' function to step-sisely reconcile base forecasts
#' 
#' @param basef basef forecasts
#' @param real dhts
#' @param lambda penalty
#' @param optimized if optimization
#' @return sw_res object holding reconciliation matrix of each step
#' @export
#' @tag single-level
stepwise_reconcile.train <- function(basef, real, lambda=0, optimized=TRUE){
  stopifnot(is(real, "dhts"))
  n = NROW(real$meta$s_mat)
  m = NCOL(real$meta$s_mat)
  domain = real$meta$domain_bts 
  if (n - m != 1){
    stop("Only two-level hierarchy is supported!")
  }
  # list storing all A_hat
  As <- list()
  bf_total <- basef[[1]]
  
  tmp <- function(x){ rbind(rep(1,x), diag(x))}
  
  for (i in 1:(m-1)){
    bf_left <- basef[[i+1]]
    if (i == m-1){
      rdomain <- domain[, (i+1):m, drop=FALSE]
      bf_right <- basef[[n]]
    } else {
      rmeta <- construct_meta(real$meta$domain_bts[,(i+1):m], tmp(m-i))
      bf_right <- marginal2Sum(basef[(i+1):n], rmeta, which = 1)
    }
    hierarchy_domain <- cbind(domain[,i], rowSums(domain[,(i+1):m, drop=FALSE]))
    # input for reconciliation of this step
    real_i <- dhts(unname(cbind(real$bts[,i], rowSums(real$bts[,(i+1):m,drop=FALSE]))), 
                   tmp(2),
                   hierarchy_domain,
                   node_names = paste0('s', c(i, i+1, paste0(i+2, '-', n))))
    basef_i <- marginal2Joint(list(bf_total, bf_left, bf_right), real_i$meta, method='ind')
    modeli <- allreconcile.train(basef_i, real_i, lambda=lambda, optimized = optimized)
    Ahat <- list(A=modeli,
                 meta=real_i$meta,
                 meta_right=rmeta)
    As[[i]] <- Ahat
    bf_total <- Joint2Marginal(reconcile(modeli, list(bf_total, bf_left, bf_right), real_i$meta), real_i$meta, 3)
  }
  structure(As, class = 'sw_res')
}

#' generic function to reconcile base forecasts
#'
#' @param x trained reconciliation model
#' @export
reconcile <- function(x, ...){
  UseMethod("reconcile", x)
}

#' function to warn irregular object
#' @param x output of reconcile.train
#' @export
reconcile.default <- function(x, ...){
  warning(paste("reconcile does not know how to handle object of class ", 
                class(x)))
}

#' step wise reconciliation 
#' 
#' @param x sw_res object
#' @param basef listing containing base probabilistic forecasts for each series
#' @export
reconcile.sw_res <- function(x, basef){
  stopifnot(is(x, 'sw_res'))
  n <- length(x)
  bf_total <- NULL
  ans <- list()
  for (i in 1:n){
    if (is.null(bf_total)) {
      bf_total <- basef[[i]]
    } else {
      bf_total <- Joint2Marginal(reconciledi, meta, 3) 
    }
    meta <- x[[i]]$meta
    rmeta <- x[[i]]$meta_right
    bf_left <- basef[[i+1]]
    if (i == n){
      bf_right <- basef[[length(basef)]]
    }else{
      bf_right <- marginal2Sum(basef[(i+1):length(basef)], rmeta, which = 1)
    }
    
    reconciledi <- reconcile(x[[i]]$A, list(bf_total, bf_left, bf_right), meta)
    ans[[i]] <- list(dist=reconciledi, meta=meta)
  }
  structure(Step2Joint(ans)$dist, class=c("coherent", "jdist", "rec"))
}


#' regular reconciliation 
#' 
#' @param x rec_mat object
#' @param basef listing containing base probabilistic forecasts for each series
#' @param domain domain of 
#' @return reconciled 
#' @export
reconcile.rec_mat <- function(x, basef, meta){
  basef <- marginal2Joint(basef, meta, method="ind")
  structure(t(x %*% t(basef)), class=c("coherent", "jdist", "rec"), meta=meta)
}

#' function to train the naive top-down method based on concurrence.
#' 
#' @param real dhts object
#' @return topdown object 
#' @export
#' @import dplyr
#' @import Matrix
topdown.train <- function(real){
  stopifnot(is(real, "dhts"))
  coherent_domain <- data.frame(unclass(real$meta$coherent_domain))
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
  
  structure(A, class = "topdown")
}

#' topdown reconcilation function
#' 
#' @param x topdown object
#' @param basef base forecasts of total level
#' @export
reconcile.topdown <- function(x, basef){
  total_f <- basef[[1]]
  structure(t(x %*% t(total_f)), class=c('coherent', 'jdist', 'td'))
}

#' function to stepwise update the joint distribution and joint domain.
#'
#' @param ans list containing joint distribution of all sub trees.
#' @return reconciled joint distribution  
#' @import dplyr
Step2Joint <- function(ans){
  
  jdist <- ans[[1]]$dist
  jmeta <- ans[[1]]$meta
  for (i in 1:(length(ans) - 1)){
    # fl: reconciled forecasts in the left hierarchy (as right node)
    fl <- Joint2Marginal(jdist, jmeta, NROW(jmeta$s_mat))
    # fr: reconciled forecasts in the right hierarchy (as total)
    rmeta <- ans[[i+1]]$meta
    fr <- Joint2Marginal(ans[[i+1]]$dist, rmeta, 1)
    adj_dist <- (fl + fr)/2
    joint_l <- coherentadj(jdist, adj_dist, jmeta, NROW(jmeta$s_mat))
    joint_r <- coherentadj(ans[[i+1]]$dist, adj_dist, rmeta, 1)
    jdist <- update_jointDist(joint_l, joint_r, jmeta, rmeta)
    jmeta <- update_meta(jmeta, rmeta)
  }
  list(meta=jmeta, dist=unname(jdist))
}

#' function to adjust joint distribution given new marginal distribution
#' 
#' @param x original joint distribution
#' @param target new marginal distribution
#' @param meta metadata of the hierarchy
#' @param which integer indicating which column in the domain 
#' the marginal distribution refers to.
#' @return adjusted joint distribution
#' @import dplyr
coherentadj <- function(x, target, meta, which){
  
  domain <- meta$coherent_domain
  time_window <- NROW(x)
  new_x <- NULL
  for (i in 1:time_window){
    tmp <- data.frame(domain=domain[,which], prob=x[i,]) %>%
      group_by(domain) %>% 
      mutate(prob = prob * target[i, as.character(pull(cur_group(), domain))] / sum(prob) ) %>%
      pull(prob)
    new_x <- rbind(new_x, tmp)
  }
  structure(unname(new_x), class = class(x))
}


addProbJitter <- function(x){
  add_x <- matrix(runif(length(x), -min(abs(x))/100, min(abs(x))/100), NROW(x), NCOL(x))
  t(apply((x + add_x), 1, function(x){abs(x)/sum(abs(x))}))
}



#' utility functions
#' @tag multiple-levels
getVIndexs <- function(distance){
  r <- dim(distance)[1]
  q <- dim(distance)[2]
  
  tmp <- (distance == matrix(apply(distance, 2, min), r, q, byrow=TRUE))
  apply(tmp, 1, which)
}

# function to construct Q matrix used in optimization
# r: cardinality of coherent domain
# pi_hat: T*q base probabilistic forecasts
construct_Q <- function(r, basef, indicators=NULL) {
  if (Matrix::rankMatrix(basef) != min(dim(basef))){
    basef <- addProbJitter(basef)
  }
  Q <- t(basef) %*% basef
  if (!corpcor::is.positive.definite(Q)){
    Q <- Q + diag(rep(1e-8, dim(Q)[1]))
  }
  if(is.null(indicators)){
    Matrix::bdiag(replicate(r, Q, simplify = FALSE))
  }else {
    Qs <- lapply(indicators, function(indicator){
      Q[indicator,][,indicator]
    })
    bdiag(Qs)
  }
}

#' utility functions
construct_D <- function(jdist, indicators, y){
  r <- NCOL(y)
  if (is.null(indicators)){
    tmp <- do.call(c, lapply(1:r, function(i){
      matrix(colSums(jdist[y[,i]==1,,drop=FALSE]), nrow = 1)
    }))
  } else {
    tmp <- do.call(c, lapply(1:r, function(i){
      matrix(colSums(jdist[y[,i]==1,,drop=FALSE])[indicators[[i]]], nrow = 1)
    }))
  }
  tmp
}

#' utility functions
construct_E <- function(distance, indicators){
  tmp1 <- do.call(c, indicators)
  tmp2 <- seq_along(tmp1)
  
  A1 <- matrix(0, NCOL(distance), length(tmp1))
  for (col in unique(tmp1)){
    A1[col, tmp2[which(tmp1 == col)]] = 1
  }
  
  A2 <- diag(length(tmp1))
  A3 <- -diag(length(tmp1))
  
  b1 <- rep(1, length(unique(tmp1)))
  b2 <- rep(0, length(tmp1))
  b3 <- rep(-1, length(tmp1))
  
  list(E=rbind(A1, A2, A3), b0=c(b1,b2,b3))
}

#' utility functions
construct_res <- function(distance, solution, indicators){
  
  res <- matrix(0, NROW(distance), NCOL(distance))
  start <- 0
  for (i in 1:length(indicators)){
    res[i, indicators[[i]]] <- solution[(start+1):(length(indicators[[i]]) + start)]
    start <- start + length(indicators[[i]])
  }
  res
}

# function to update joint distribution
update_jointDist <- function(ldist, rdist, lmeta, rmeta){
  ld <- data.frame(unclass(lmeta$coherent_domain))
  rd <- data.frame(unclass(rmeta$coherent_domain))
  colnames(ld)[dim(ld)[2]] <- "key"
  colnames(rd)[1] <- "key" 
  time_window <- dim(ldist)[1]
  new_dist <- NULL
  for (i in 1:time_window){
    ld$lprob <- ldist[i,]
    new_dist <- mutate(rd, prob=rdist[i,]) %>% group_by(key) %>%
      mutate(ratio=prob/sum(prob)) %>%
      right_join(ld, by="key") %>%
      mutate(prob = lprob * ratio) %>%
      pull(prob) %>%
      rbind(new_dist, .)
  }
  structure(as.matrix(unname(new_dist)), class=class(ldist))
}

update_meta <- function(l, r){
  
  new_domain <- cbind(l$domain_bts[,1:(NCOL(l$domain_bts) - 1)], r$domain_bts)
  construct_meta(new_domain,
                 rbind(rep(1, NCOL(new_domain)), diag(NCOL(new_domain))),
                 node_names = c(colnames(l$coherent_domain)[1:(NCOL(l$coherent_domain) - 1)],
                                colnames(r$coherent_domain)[2:NCOL(r$coherent_domain)]))
  
  # l <- data.frame(unclass(l$coherent_domain))
  # r <- data.frame(unclass(r$coherent_domain))
  # colnames(l)[dim(l)[2]] <- "key"
  # colnames(r)[1] <- "key"
  # ans <- right_join(r, l, by="key")
  # ans <- ans %>% select(-c("key")) %>% 
  #   select(c(3:dim(.)[2], 1:2))
  # colnames(ans) <- 1:dim(ans)[2]
  
  # l$meta$coherent_domain <- unname(structure(as.matrix(ans), class="coherent_domain"))
}

