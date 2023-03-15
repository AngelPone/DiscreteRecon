#' reconcile base probabilistic forecasts
#' 
#' Method for reconciling base probabilistic forecasts 
#' @param bf_train incoherent joint base probabilistic forecasts
#' @param obs_train dhts object holding observed series
#' @param hier dhier object
#' @param lambda penalty terms of cost, default is 0. not working if optimized=
#' TRUE
#' @param optimized indicating whether optimized. i.e., 
#' only moving probabilities to nearest coherent points.
#' @return transformation matrix
fr_train <- function(hier,
                     bf_train, 
                     obs_train,
                     lambda=0,
                     optimized=TRUE) {
  
  bf_train <- marginal2Joint(bf_train, hier)
  real_dummy <- cons_realDummy(obs_train, hier)
  
  cdomain <- hier$coherent_domain
  idomain <- hier$incoherent_domain
  
  r <- NROW(cdomain)
  q <- NROW(idomain)
  time_window <- NROW(real_dummy)
  costs <- cal_costeMatrix(hier)
  
  if (optimized){
    indicators <- getVIndexs(costs)
    Dmat <- construct_Q(r, bf_train, indicators)
    vec_c <- do.call(c, lapply(1:r, function(i){costs[i,indicators[[i]]]}))
    dvec <- - 2 * construct_D(bf_train, indicators, real_dummy)
    n <- sum(sapply(indicators, length))
    A1 <- construct_E(costs, indicators)
  } else {
    n <- r * q
    A1 <- do.call(cbind, replicate(r, Diagonal(q)))
    Dmat <- construct_Q(r, bf_train)  
    vec_c <- as.vector(t(costs))
    dvec <- time_window * lambda * vec_c - 2 * construct_D(bf_train, NULL, real_dummy)
  }
  vs <- Variable(n)
  constraints <- list(A1 %*% vs == 1, vs >= 0, vs <= 1)
  objective <- Minimize(t(dvec) %*% vs + quad_form(vs, Dmat))
  problem <- Problem(objective, constraints)
  solution <- solve(problem)
  if (optimized){
    Ahat <- construct_res(costs, solution$getValue(vs), indicators)
  }else{
    Ahat <- matrix(solution$getValue(vs), r, q, byrow = TRUE)
  }
  apply(Ahat, 2, function(x){x/sum(x)})
}

#' validator for input observations
#' 
#' @param x time series
#' @param hier hierarchy
validate_obs <- function(x, hier){
  stopifnot(NCOL(x) == NCOL(hier$s_mat))
  if (any(apply(x, 2, max) > hier$domain[2,],
          apply(x, 2, min) > hier$domain[1,])) {
    stop("The observations should be in the range of the domain.")
  }
}

#' validator for input base forecasts
#' 
#' @param x base probabilistic forecasts
#' @param hier hierarchy
validate_bf <- function(x, hier){
  stopifnot(is.list(x),
            length(x) == hier$quantities['n'])
  ranges <- hier$s_mat %*% t(hier$domain)
  
  for (i in 1:NROW(ranges)){
    f <- x[[i]]
    if (any(colnames(f) != ranges[i, 1] : ranges[i, 2])) {
      stop(paste0("colnames of the ", i, "th series should equal to its domain"))
    }
  }
}

#' function to train discrete reconciliation models
#' 
#' @param hier `dhier` object representing the hierarchy.
#' @param method Available methods include "dfr", "sdfr", "bu", "td".
#' @param bf_train list containing base probabilistic forecasts of each series.
#' It should be provided for method "dfr" and "sdfr".
#' @param obs_train Past observations used for train the model.
#' It should be provided for all methods except "bu".
#' @param lambda penalty value used in "dfr"
#' @param optimized boolean value used in "dfr". If TRUE, the elements of 
#' non-nearest points in the reconciliation matrix will be zero.
#' @return the reconciliation model
#' @export
dfr <- function(hier,
                method=c("bu", "td", "dfr", "sdfr"),
                obs_train=NULL, 
                bf_train=NULL,
                lambda=0, 
                optimized=TRUE){
  stopifnot(is.dhier(hier))
  method <- match.arg(method)
  
  if (!is.null(obs_train)) validate_obs(obs_train, hier)
  if (!is.null(bf_train)) validate_bf(bf_train, hier)
  
  output <- switch(method,
                   bu = bu_train(hier),
                   td = td_train(hier, obs_train),
                   dfr = fr_train(hier, bf_train, obs_train, lambda, optimized),
                   sdfr = sfr_train(hier, bf_train, obs_train, lambda, optimized))
  structure(list(A = output, hier=hier, method=method), class = "dfr")
}

#' function to stepsise reconcile base forecasts
#' 
#' @param hier dhier object
#' @param bf_train list containing base probabilistic forecasts of each series.
#' It should be provided for method "dfr" and "sdfr".
#' @param obs_train Past observations used for train the model
#' @param lambda penalty value
#' @param optimized boolean
#' @return reconciliation matrices
sfr_train <- function(hier, bf_train, obs_train, lambda=0, optimized=TRUE){
  n = NROW(hier$s_mat)
  m = NCOL(hier$s_mat)
  if (m <= 2) {
    stop("SDFR can only be used for hierarchies with more than two bottomlevel series.")
  }
  domain = hier$domain
  if (n - m != 1){
    stop("Only two-level hierarchy is supported!")
  }
  # list storing all A_hat
  As <- list()
  bf_total <- bf_train[[1]]
  
  tmp <- function(x){ rbind(rep(1,x), diag(x)) }
  
  for (i in 1:(m-1)){
    bf_left <- bf_train[[i+1]]
    if (i == m-1){
      rdomain <- domain[, (i+1):m, drop=FALSE]
      bf_right <- bf_train[[n]]
    } else {
      rmeta <- dhier(tmp(m-i), hier$domain[,(i+1):m])
      bf_right <- marginal2Sum(bf_train[(i+1):n], rmeta, which = 1)
    }
    hierarchy_domain <- cbind(domain[,i], rowSums(domain[,(i+1):m, drop=FALSE]))
    # input for reconciliation of this step
    hier_i <- dhier(tmp(2),
                    hierarchy_domain,
                    node_names = paste0('s', c(i, i+1, paste0(i+2, '-', n))))
    modeli <- dfr(hier_i, 
                  method = "dfr",
                  bf_train = list(bf_total, bf_left, bf_right), 
                  obs_train = unname(cbind(obs_train[,i], rowSums(obs_train[,(i+1):m,drop=FALSE]))), 
                  lambda=lambda, optimized = optimized)
    Ahat <- list(model=modeli,
                 meta_right=rmeta)
    As[[i]] <- Ahat
    bf_total <- Joint2Marginal(reconcile(modeli, list(bf_total, bf_left, bf_right)), hier_i, 3)
  }
  As
}

#' print function for dfr object
#' 
#' @param x dfr object
#' @export
print.dfr <- function(x){
  cat(paste0("trained ", 
             switch(x$method,
                    td = "Top-Down",
                    bu = "Bottom-Up",
                    dfr = "DFR",
                    sdfr = "SDFR"),
             " reconciliation model for \n"))
  print.dhier(x$hier)
  invisible(x)
}

#' generic function to reconcile base forecasts
#'
#' @param mdl trained reconciliation model
#' @param new_bf base forecasts
#' @return reconciled forecasts
#' @export
reconcile <- function(mdl, new_bf){
  stopifnot("dfr" %in% class(mdl))
  validate_bf(new_bf, mdl$hier)
  if (is.list(mdl$A)){
    reconcile_sdfr(mdl, new_bf)
  } else {
    reconcile_dfr(mdl, new_bf)
  }
}

#' stepwise reconciliation 
#' 
#' @param mdl trained stopwise reconciliation model
#' @param new_bf base forecasts
#' @return reconciled forecasts
reconcile_sdfr <- function(mdl, new_bf){
  n <- length(mdl$A)
  bf_total <- NULL
  ans <- list()
  for (i in 1:n){
    if (is.null(bf_total)) {
      bf_total <- new_bf[[i]]
    } else {
      bf_total <- Joint2Marginal(reconciledi, meta, 3) 
    }
    meta <- mdl$A[[i]]$model$hier
    rmeta <- mdl$A[[i]]$meta_right
    bf_left <- new_bf[[i+1]]
    if (i == n){
      bf_right <- new_bf[[length(new_bf)]]
    }else{
      bf_right <- marginal2Sum(new_bf[(i+1):length(new_bf)], rmeta, which = 1)
    }
    
    reconciledi <- reconcile(mdl$A[[i]]$model, list(bf_total, bf_left, bf_right))
    ans[[i]] <- list(dist=reconciledi, hier=meta)
  }
  Step2Joint(ans)
}


#' regular reconciliation 
#' 
#' @param mdl trained reconciliation model
#' @param bf base forecasts
#' @return reconciled forecasts
reconcile_dfr <- function(mdl, bf){
  bf <- marginal2Joint(bf, mdl$hier)
  t(as.matrix(mdl$A) %*% t(bf))
}

#' function to stepwise update the joint distribution and joint domain.
#'
#' @param ans list containing joint distribution of all sub trees.
#' @return reconciled joint distribution  
Step2Joint <- function(ans){
  
  jdist <- ans[[1]]$dist
  jmeta <- ans[[1]]$hier
  for (i in 1:(length(ans) - 1)){
    # fl: reconciled forecasts in the left hierarchy (as right node)
    fl <- Joint2Marginal(jdist, jmeta, NROW(jmeta$s_mat))
    # fr: reconciled forecasts in the right hierarchy (as total)
    rmeta <- ans[[i+1]]$hier
    fr <- Joint2Marginal(ans[[i+1]]$dist, rmeta, 1)
    adj_dist <- (fl + fr)/2
    joint_l <- coherentadj(jdist, adj_dist, jmeta, NROW(jmeta$s_mat))
    joint_r <- coherentadj(ans[[i+1]]$dist, adj_dist, rmeta, 1)
    jdist <- update_jointDist(joint_l, joint_r, jmeta, rmeta)
    jmeta <- update_meta(jmeta, rmeta)
  }
  unname(jdist)
}

#' function to adjust joint distribution given new marginal distribution
#' 
#' @param x original joint distribution
#' @param target new marginal distribution
#' @param hier metadata of the hierarchy
#' @param which integer indicating which column in the domain 
#' the marginal distribution refers to.
#' @return adjusted joint distribution
coherentadj <- function(x, target, hier, which){
  domain <- hier$coherent_domain
  for (j in unique(domain[,which])){
    idx <- domain[,which] == j
    x[,idx] <- x[,idx] * target[,paste0(j)] / rowSums(x[,idx, drop=FALSE]) 
  }
  x
}


#' utility functions
getVIndexs <- function(distance){
  r <- dim(distance)[1]
  q <- dim(distance)[2]
  
  tmp <- (distance == matrix(apply(distance, 2, min), r, q, byrow=TRUE))
  apply(tmp, 1, which)
}

#' function to construct Q matrix used in optimization
#' 
#' @param r cardinality of coherent domain
#' @param basef base forecasts 
#' @param indicators indicators
#' @return Q matrix
construct_Q <- function(r, basef, indicators=NULL) {
  Q <- t(basef) %*% basef
  if (!is.positive.definite(Q)){
    Q <- Q + diag(rep(1e-8, dim(Q)[1]))
  }
  if(is.null(indicators)){
    bdiag(replicate(r, Q, simplify = FALSE))
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
      colSums(jdist[y[,i]==1,,drop=FALSE])
    }))
  } else {
    tmp <- do.call(c, lapply(1:r, function(i){
      colSums(jdist[y[,i]==1,,drop=FALSE])[indicators[[i]]]
    }))
  }
  tmp
}

#' utility functions
construct_E <- function(distance, indicators){
  tmp1 <- do.call(c, indicators)
  tmp2 <- seq_along(tmp1)
  
  A1 <- Matrix::spMatrix(NCOL(distance), length(tmp1))
  for (col in unique(tmp1)){
    A1[col, tmp2[which(tmp1 == col)]] = 1
  }
  as(A1, "CsparseMatrix")
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

#' function to update joint distribution
update_jointDist <- function(ldist, rdist, lmeta, rmeta){
  new_hier <- update_meta(lmeta, rmeta)
  domain <- new_hier$coherent_domain
  key <- rowSums(domain[,(NCOL(domain)-1) : NCOL(domain)])
  ldomain <- lmeta$coherent_domain
  rdomain <- rmeta$coherent_domain
  output <- matrix(0, NROW(ldist), NROW(domain))
  for (j in unique(key)){
    idx <- which(key == j)
    lidx <- which(ldomain[,NCOL(ldomain)] == j)
    ridx <- which(rdomain[,1] == j)
    if (length(ridx) == 1){
      output[,idx] <- ldist[,lidx]
    } else {
      ratios <- rdist[,ridx,drop=FALSE] / rowSums(rdist[,ridx, drop=FALSE])
      for (curd in idx){
        lidx_i <- lidx[get_equal_idx(ldomain[lidx, 1:(NCOL(ldomain) - 1)], 
                                  domain[curd, 1:(NCOL(domain)-2)])]
        ridx_i <- get_equal_idx(rdomain[ridx, 2:3], 
                                domain[curd, (NCOL(domain)-1):NCOL(domain)])
        output[,curd] <- ldist[,lidx_i] * ratios[,ridx_i]
      }
    }
  }
  output
}

#' utility
get_equal_idx <- function(ldf, row){
  for (i in 1:NROW(ldf)){
    if(all(ldf[i,] == row)) break
  }
  i
}

#' utility function
update_meta <- function(l, r){
  
  new_domain <- cbind(l$domain[,1:(NCOL(l$domain_bts) - 1)], r$domain)
  dhier(rbind(rep(1, NCOL(new_domain)), diag(NCOL(new_domain))),
        new_domain,
        node_names = c(colnames(l$coherent_domain)[1:(NCOL(l$coherent_domain) - 1)],
                       colnames(r$coherent_domain)[2:NCOL(r$coherent_domain)]))
}



#' function to train bottom-up model
#' 
#' @param x dhier object
#' @return reconciliation matrix
bu_train <- function(x){
  d <- x$quantities['q'] / x$quantities['r']
  Matrix::bdiag(rep(list(matrix(1, 1, d)), x$quantities['r']))
}


#' function to train top-down model
#' 
#' @param x dhier object
#' @param obs_train observations used for train
#' @return reconciliation matrix
td_train <- function(x, obs_train){
  n <- x$quantities['n']
  m <- x$quantities['m']
  r <- x$quantities['r']
  q <- x$quantities['q']

  concurrence <- vector("numeric", r)
  for (i in 1:NROW(obs_train)){
    for (j in 1:r){
      if (all(obs_train[i,] == x$coherent_domain[j, (n-m+1):n])){
        concurrence[j] = concurrence[j] + 1
        break
      }
    }
  }
  
  probs <- unsplit(lapply(split(concurrence, x$coherent_domain[,1]), 
                          function(x){ x/sum(x) }),
                   x$coherent_domain[,1])
  output <- Matrix::spMatrix(r, q)
  for (j in 1:r){
    output[j, which(x$incoherent_domain[,1] == x$coherent_domain[j, 1])] = 
      probs[j]
  }
  output
}

