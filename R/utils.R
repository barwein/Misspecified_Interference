
# library(igraph)
# library(data.table)
# library(mgcv)
# library(parallel)
# library(gtools)

# Basic R functions -----------------------------------------------------------

#' Compute exposure values in the four level thresholds model
#' @param A adj. matrix of a network
#' @param Z treatment vector
#' @param threshold vector of unit-level thresholds
#'
#' @returns vector of observed exposures
#' @export
generate_exposures_threshold <- function(A, Z, threshold,deg=NULL){
  # Compute num neighbors
  if(is.null(deg)){
    deg <- apply(A, 1, sum)
  }
  # Compute signs
  indirect_sign <- as.numeric((A %*% Z) > threshold*deg)

  # Save exposures
  exposures <- vector(length = nrow(A))

  exposures[Z*indirect_sign==1] <- "c11"
  exposures[(1-Z)*indirect_sign==1] <- "c01"
  exposures[Z*(1-indirect_sign)==1] <- "c10"
  exposures[(1-Z)*(1-indirect_sign)==1] <- "c00"

  return(exposures)
}

#' Return a matrix of exposures
#' @noRd
exposures_under_network_list <- function(Z, A.list, exposure_func, param.list){
  # for m adj. matrices in A.list and n units, this function compute the nXm matrix of exposures under
  # each adj. matrix for treatment vector Z
    return(
      # sapply(A.list, function(A){generate_exposures_threshold(A=A,Z=Z,threshold=threshold)})
      sapply(A.list, function(A){do.call(exposure_func,list(A,Z,unlist(param.list)))})
      )}

#' check which units have exposure ck under all networks
#' @noRd
specific_exposure_under_network_list <- function(exposures_mat, ck){
  # for some exposures matrix, this function return a binary (1/0) vector
  # of units that have exposure ck under all adj. matrics (rows in exposure_mat)
  return(
    apply(exposures_mat,1,function(x){as.numeric(all(x==ck))})
         )}

#' A function that returns a list of nXr indicator matrices (one for each exposures),
#' that indicate if unit i in sample r is exposed to some ck under all A \in A.list
#' The output is used for computation of exposures probabilities
#' @noRd
compute_prob_indicator <- function(R, n, Pz_function, pz_func_args, A.list, exposures_vec, exposure_func, exposure_func_args){
  # Init a list of nXR indicator matrices
  Ind.mat.list <- vector("list",length(exposures_vec))
  for (i in seq(length(exposures_vec))) {
    assign(paste0("I.",exposures_vec[i]), matrix(NA, nrow = n, ncol = R))
    Ind.mat.list[[i]] <- get(paste0("I.",exposures_vec[i]))
  }
  names(Ind.mat.list) <- exposures_vec

  # Sample Z ~ PZ R times and update the indicators matrices
  for (r in seq(R)) {
    # cur.Z <- Pz_function(n=n, p=pz)
    cur.Z <- do.call(Pz_function, pz_func_args)
    # cur.exposure.mat <- exposures_under_network_list(Z=cur.Z,
    #                                                  A.list=A.list,
    #                                                  threshold=threshold)
    cur.exposure.mat <- exposures_under_network_list(Z=cur.Z,
                                                     A.list=A.list,
                                                     exposure_func = exposure_func,
                                                     param.list = exposure_func_args)
    for (ck in exposures_vec) {
      # Update indicשtor matrices for each exposures and r
      Ind.mat.list[[ck]][,r] <- specific_exposure_under_network_list(exposures_mat=cur.exposure.mat,
                                                                     ck=ck)
    }
  }
  return(Ind.mat.list)
}

#' From the list of indicator matrices, vector of exposures and relevant contrast,
#' the function return a list of all relevant probability matrices
#' @noRd
compute_prob_matrices <- function(ind.mat.list, exposures_contrast, exposures_vec){
  n <- nrow(ind.mat.list[[1]])
  R <- ncol(ind.mat.list[[1]])
  # Init list that will contain all relevant probability matrices
  P_list <- vector(mode = "list", length = length(exposures_vec) + length(exposures_contrast))

  # Update P.k matrices
  for (i in seq(length(exposures_vec))) {
    ind.mat <- ind.mat.list[[exposures_vec[i]]]
    # Use additive smoothing for the diagonal values
    P_list[[i]] <- (ind.mat %*% t(ind.mat) + diag(1,n)) / (R+1)
  }

  if(length(exposures_contrast) > 0){
  # Update P.kl matrices
    for (j in seq(length(exposures_contrast))) {
      ck <- exposures_contrast[[j]][1]
      cl <- exposures_contrast[[j]][2]
      ind.mat.k <- ind.mat.list[[ck]]
      ind.mat.l <- ind.mat.list[[cl]]
      P_list[[length(exposures_vec)+j]] <- (ind.mat.k %*% t(ind.mat.l)) / R
    }
  }
  names(P_list) = c(exposures_vec,
                   sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])}))

  return(P_list)
}


#' Wrapper function that take as input all the adj. matrix list, R (number of iteration),
#' exposures, exposures contrasts list, and Z sampler,
#' and compute the list of probability matrices for the later use by MR estimator
#' @noRd
Get_prob_matrices_list <- function(R,
                                   n,
                                   Pz_function,
                                   pz_func_args,
                                   A.list,
                                   exposures_contrast, exposures_vec,
                                   exposure_func,
                                   exposure_func_args,
                                   # threshold,
                                   Palluck.et.al = FALSE,
                                   schid = NULL){

  # First, get the indicator matrix list (list with length equal to exposures_vec,
  #                                         one matrix for each exposure)
  if(Palluck.et.al){
    ind.mat.list <- compute_prob_indicator_palluck(R = R, n = n,
                                                   Pz_function = Pz_function,
                                                   pz_func_args = pz_func_args,
                                                   # schid = schid,
                                                   A.list = A.list,
                                                   exposures_vec = exposures_vec)
  } else{
    ind.mat.list <- compute_prob_indicator(R = R, n = n,
                                           Pz_function = Pz_function,
                                           pz_func_args = pz_func_args,
                                           A.list = A.list,
                                           exposures_vec = exposures_vec,
                                           exposure_func = exposure_func,
                                           # threshold = threshold)
                                           exposure_func_args = exposure_func_args)
  }

  # Compute the relevant prob. matrices
  Prob_matrices_list <- compute_prob_matrices(ind.mat.list = ind.mat.list,
                                              exposures_contrast = exposures_contrast,
                                              exposures_vec = exposures_vec)

  return(Prob_matrices_list)
}

#' Function that estimate the mean potential outcome using the MR estimator
#' @noRd
compute_MR_PO_estimator <- function(Y.obs, expos.obs.mat, ck, P.k){

  n <- length(Y.obs)
  # Get indicator of exposure ck vector
  ind.k <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = ck)

  # Get probabilities from the diagonal of P.k
  Pii <- diag(P.k)

  # Estimate using both HT and Hajek

  HT_esti <- (n^{-1})*(sum(ind.k*Y.obs/Pii))

  hajek_esti <- (sum(ind.k*Y.obs/Pii)) / (sum(ind.k/Pii))

  return(list(ht_esti = HT_esti, hajek_esti = hajek_esti))
}

#' Function that compute the variance of the mean potential outcome  MR estimator
#' @noRd
compute_var_MR_PO_estimator <- function(Y.obs, expos.obs.mat, ck, P.k, estimate.n = FALSE){

  n <- length(Y.obs)
  # Get indicator of exposure ck vector
  ind.k <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = ck)

  # Get probabilities from the diagonal of P.k
  Pii <- diag(P.k)

  # Whether to estimate n or use given n (HT vs Hajek)
  n.hat <- ifelse(estimate.n,sum(ind.k/Pii),n)

  # First sum of the formula
  sum.1 <- (n.hat^{-2})*(sum(ind.k*
                        (1-Pii)*
                        ((Y.obs/Pii)^2)
                        ))

  sum.2 = 0; sum.3 = 0

  # Iteratively update second and third of the formula
  for (i in seq(n)) {
    for (j in setdiff(seq(n),i)){
      if(P.k[i,j] != 0 & ind.k[i]*ind.k[j]==1){
      # if(ind.k[i]*ind.k[j]==1){
        cur.sum = ((P.k[i,j] - Pii[i]*Pii[j])/(P.k[i,j]))*
                   ((Y.obs[i]*Y.obs[j])/(Pii[i]*Pii[j]))
        sum.2 = sum.2 + cur.sum
      }
      if(P.k[i,j] == 0){
        cur.sum = (ind.k[i]*(Y.obs[i]^2)/(2*Pii[i])) +
                  (ind.k[j]*(Y.obs[j]^2)/(2*Pii[j]))
        sum.3 = sum.3 + cur.sum
      }
    }
  }
  # Scale the sums
  sum.2 = (n.hat^{-2})*sum.2
  sum.3 = (n.hat^{-2})*sum.3
  # Return variance estimator
  return(sum.1 + sum.2 + sum.3)
}


#' Compute the COV term of NMR
#' @noRd
compute_cov_MR_PO_estimator <- function(Y.obs, expos.obs.mat, ck, cl,
                                        P.k, P.l, P.kl,
                                        estimate.n = FALSE){

  n <- length(Y.obs)
  # Get indicator of exposure ck vector
  ind.k <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = ck)
  ind.l <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = cl)

  # Whether to estimate n or use given n (HT vs Hajek)
  # n.hat <- ifelse(estimate.n,sum(ind.k/Pii),n)
  n.hat <- n

  sum.1 = 0; sum.2 = 0

  for (i in seq(n)) {
    for (j in seq(n)) {
      if(P.kl[i,j] > 0 & ind.k[i]*ind.l[j] == 1 & i != j){
        cur.sum = ((P.kl[i,j] - P.k[i,i]*P.l[j,j])/(P.kl[i,j]))*
                  ((Y.obs[i]*Y.obs[j])/(P.k[i,i]*P.l[j,j]))
        sum.1 = sum.1 + cur.sum
      }
      if(P.kl[i,j] == 0){
        cur.sum = (ind.k[i]*(Y.obs[i]^2)/(2*P.k[i,i])) +
                  (ind.l[j]*(Y.obs[j]^2)/(2*P.l[j,j]))
        sum.2 = sum.2 + cur.sum
      }
    }
  }
  sum.1 = (n.hat^{-2})*sum.1
  sum.2 = (n.hat^{-2})*sum.2

  return(sum.1 - sum.2)
}


#' Get linear residualized values from Hajek estimator
#' only for units with exposure ck under all adj. matrix list.
#' @noRd
Hajek_residualized_Y <- function(Y.obs, ck, cl, expos.obs.mat, hajek_esti.ck, hajek_esti.cl){

  ind.k <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = ck)
  ind.l <- specific_exposure_under_network_list(exposures_mat = expos.obs.mat, ck = cl)
  y.resid <- ind.k*(Y.obs - hajek_esti.ck) + ind.l*(Y.obs - hajek_esti.cl)
  return(y.resid)
}


#' Wrapper function that compute point estimate and variance of causal effect \tau(ck,cl)
#' @noRd
MR_CE_wrapper <- function(Y.obs,
                          expos.obs.mat,
                          ck,
                          cl,
                          P.k,
                          P.l,
                          P.kl,
                          estimate.n = FALSE,
                          compute.variance = TRUE){
  # Get point estimates using HT and Hajek
  mu.ck <- compute_MR_PO_estimator(Y.obs = Y.obs,
                                   expos.obs.mat = expos.obs.mat,
                                   ck = ck,
                                   P.k = P.k)

  mu.cl <- compute_MR_PO_estimator(Y.obs = Y.obs,
                                   expos.obs.mat = expos.obs.mat,
                                   ck = cl,
                                   P.k = P.l)

  ht_ce <- mu.ck$ht_esti - mu.cl$ht_esti
  hajek_ce <- mu.ck$hajek_esti - mu.cl$hajek_esti

  if(compute.variance){
    # Get variance estimator for HT
    var.ht.ck <- compute_var_MR_PO_estimator(Y.obs = Y.obs,
                                              expos.obs.mat = expos.obs.mat,
                                              ck = ck,
                                              P.k = P.k,
                                             estimate.n = estimate.n)

    var.ht.cl <- compute_var_MR_PO_estimator(Y.obs = Y.obs,
                                              expos.obs.mat = expos.obs.mat,
                                              ck = cl,
                                              P.k = P.l,
                                             estimate.n = estimate.n)

    cov.ht.ckcl <- compute_cov_MR_PO_estimator(Y.obs = Y.obs,
                                                expos.obs.mat = expos.obs.mat,
                                                ck = ck,
                                                cl = cl,
                                                P.k = P.k,
                                                P.l = P.l,
                                                P.kl = P.kl,
                                               estimate.n = estimate.n)

    var_ht_ce <- var.ht.ck + var.ht.ck - 2*cov.ht.ckcl

    # Get variance estimator for Hajek
    # First need to residualize Y.obs and then use mentioned variance estimator
    Y.resid <- Hajek_residualized_Y(Y.obs = Y.obs,
                                       ck = ck,
                                       cl = cl,
                                    expos.obs.mat = expos.obs.mat,
                                       hajek_esti.ck = mu.ck$hajek_esti,
                                       hajek_esti.cl = mu.cl$hajek_esti)

    var.hajek.ck <- compute_var_MR_PO_estimator(Y.obs = Y.resid,
                                                 expos.obs.mat = expos.obs.mat,
                                                 ck = ck,
                                                 P.k = P.k,
                                                estimate.n = estimate.n)

    var.hajek.cl <- compute_var_MR_PO_estimator(Y.obs = Y.resid,
                                                 expos.obs.mat = expos.obs.mat,
                                                 ck = cl,
                                                 P.k = P.l,
                                                estimate.n = estimate.n)

    cov.hajek.ckcl <- compute_cov_MR_PO_estimator(Y.obs = Y.resid,
                                                   expos.obs.mat = expos.obs.mat,
                                                   ck = ck,
                                                   cl = cl,
                                                   P.k = P.k,
                                                   P.l = P.l,
                                                   P.kl = P.kl,
                                                  estimate.n = estimate.n)

    var_hajek_ce <- var.hajek.ck + var.hajek.cl - 2*cov.hajek.ckcl
  }
  if(!compute.variance){
    var_ht_ce = var_hajek_ce = NA
  }
  # Return results
  return(list(ht_ce = ht_ce, hajek_ce = hajek_ce,
              var_ht_ce = var_ht_ce, var_hajek_ce = var_hajek_ce))
}

#' Wrapper for CE estimation with NMR
#' @noRd
MR_CE_estimator <- function(Z.obs,
                            Y.obs,
                            A.list,
                            exposures_contrast,
                            exposures_vec,
                            Prob_matrices_list,
                            exposure_func,
                            exposure_func_args,
                            estimate.n = FALSE,
                            compute.variance = TRUE){

  # expos.obs.mat <- exposures_under_network_list(Z = Z.obs, A.list = A.list, threshold = threshold)
  expos.obs.mat <- exposures_under_network_list(Z = Z.obs, A.list = A.list,
                                                exposure_func = exposure_func, param.list = exposure_func_args)

  ce.results <- vector("list", length(exposures_contrast))
  names(ce.results) <- sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])})

  for (expos in exposures_contrast) {
    ck <- expos[1]
    cl <- expos[2]
    cntrast_name <- paste0(ck,"-",cl)
    P.k <- Prob_matrices_list[[ck]]
    P.l <- Prob_matrices_list[[cl]]
    P.kl <- Prob_matrices_list[[cntrast_name]]
    ce_estimate <- MR_CE_wrapper(Y.obs = Y.obs,
                                 expos.obs.mat = expos.obs.mat,
                                 ck = ck,
                                 cl = cl,
                                 P.k = P.k,
                                 P.l = P.l,
                                 P.kl = P.kl,
                                 estimate.n = estimate.n,
                                 compute.variance = compute.variance)
    ce.results[[cntrast_name]] <- ce_estimate
  }
  return(ce.results)
}



# PBA ----------------------------------------------------


#' Generate Q matrix of SBM
#' @noRd
generate_Q_matrix <- function(N_clusters,
                              within_prob_vec,
                              between_prob_vec){
  Q <- matrix(NA, N_clusters, N_clusters)
  diag(Q) <- within_prob_vec
  Q[lower.tri(Q)] <- between_prob_vec
  Q <- as.matrix(Matrix::forceSymmetric(Q, uplo = "L"))
  return(Q)
}

exp_spatial_distance <- function(dist_vec, theta){
  return(exp(-dist_vec/(mean(dist_vec)*theta)))
}


#' Sample treatment with Bernoulli CRT
#' @noRd
Z_ber_clusters <- function(N_clusters,
                           N_each_cluster_vec,
                           p){

  treated_cluster_indicator = rep(0, N_clusters)
  treated_cluster_idx = sample(x = seq(N_clusters),
                               size = N_clusters %/% (1/p),
                               replace = FALSE)
  treated_cluster_indicator[treated_cluster_idx] <- 1
  unit_level_treatment_vec = rep(treated_cluster_indicator, N_each_cluster_vec)
  return(unit_level_treatment_vec)
}

#' Generate PO as in the paper
#' @param exposures vector of observed exposures
#' @param base.noise baseline potential outcome
#'
#' @export
generate_po <- function(exposures, base.noise = NULL){

  if(is.null(base.noise)){
    n <- length(exposures)
    y.00 <- runif(n = n, min = 0.5, max = 1.5)
  } else{
    y.00 <- base.noise
  }
  Y <- (y.00+1)*(exposures=="c11") +
    (y.00+0.25)*(exposures=="c01") +
    (y.00+0.5)*(exposures=="c10") +
    y.00*(exposures=="c00")

  return(Y)
}

#' Create perturb betwork from A.sp via the ``edges_prob_func" selected by the user
#' @noRd
generate_perturbed_network <- function(N_units,
                                       A.sp,
                                       edges_prob_func,
                                       prior_func_list,
                                       prior_func_args_list,
                                       X.obs = NULL){
  # Sample from prior
  theta.length <- length(prior_func_list)
  theta.vec <- vector("numeric",theta.length)
  for (i in seq(theta.length)){
    theta.vec[i] <- do.call(prior_func_list[[i]],prior_func_args_list[[i]])
  }
  # Generate P_\theta prob matrix (lower tri matrix with upper tri = 0 or NA)
  if(is.null(X.obs)){
    prob.lower.tri.matrix <- do.call(edges_prob_func, list(theta.vec, A.sp))
  } else{
    prob.lower.tri.matrix <- do.call(edges_prob_func, list(theta.vec, A.sp, X.obs))
  }
  # Sample edges
  edges.prob.vec <- prob.lower.tri.matrix[lower.tri(prob.lower.tri.matrix)]
  edges.lower.tri <- rbinom(n = length(edges.prob.vec), size = 1, prob = edges.prob.vec)
  # Save as A.star
  A.star <- matrix(NA, N_units, N_units)
  A.star[lower.tri(A.star)] <- edges.lower.tri
  diag(A.star) <- 0
  # Make it symmetric (undirected)
  A.star <- Matrix::forceSymmetric(A.star, uplo = "L")
  return(A.star)
}


#' One iteration of PBA
#' @noRd
Single_PBA_general_iteration <- function(N_units,
                                         A.sp,
                                         edges_prob_func,
                                         prior_func_list,
                                         prior_func_args_list,
                                         Z.obs,
                                         Y.obs,
                                         X.obs = NULL,
                                         Pz_function,
                                         pz_func_args,
                                         exposures_vec,
                                         exposures_contrast,
                                         exposures_thresholds = NULL){

  # Generate A.star (perturbed version of A.sp)
  cur.A.star <- generate_perturbed_network(N_units = N_units,
                                           A.sp = A.sp,
                                           edges_prob_func = edges_prob_func,
                                           prior_func_list = prior_func_list,
                                           prior_func_args_list = prior_func_args_list,
                                           X.obs = X.obs)


  if(is.null(exposures_thresholds)){exposures_thresholds = rep(0, N_units)}

  # Get prob. matrices for current censored network
  prob.mat <- Get_prob_matrices_list(
    # R = 10^4,
    R = 10,
    n = N_units,
    Pz_function = Pz_function,
    pz_func_args = pz_func_args,
    A.list = list(cur.A.star),
    exposures_contrast = exposures_contrast,
    exposures_vec = exposures_vec,
    threshold = exposures_thresholds)


  cur.CE.estimates <- MR_CE_estimator(Z.obs = Z.obs,
                                      Y.obs = Y.obs,
                                      A.list = list(cur.A.star),
                                      exposures_contrast = exposures_contrast,
                                      exposures_vec = exposures_vec,
                                      Prob_matrices_list = prob.mat,
                                      exposure_func = generate_exposures_threshold,
                                      exposure_func_args = list(A=cur.A.star,
                                                                Z=Z.obs,
                                                                threshold = exposures_thresholds))

  cur.CE.estimates <- data.table(do.call(rbind,cur.CE.estimates))

  # Add random error via normal approximation
  cur.CE.estimates$ht_ce_w_re <- mapply(rnorm,
                                        n=1,
                                        mean=cur.CE.estimates$ht_ce,
                                        sd=sqrt(cur.CE.estimates$var_ht_ce))

  cur.CE.estimates$hajek_ce_w_re <- mapply(rnorm,
                                           n=1,
                                           mean=cur.CE.estimates$hajek_ce,
                                           sd=sqrt(cur.CE.estimates$var_hajek_ce))

  cur.CE.estimates$ce_contrast <- sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])})

  return(cur.CE.estimates)
}

