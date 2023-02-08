
# source("R/Aux_functions.R")


#' Sensitivity analysis (SA) for misspecification of the network interference structure
#'
#' @description
#' Return the causal estimates from the obtained by perturbing the specified network
#'
#' @details
#' Examine how sensitive causal estimates are to deviations of the specified network;
#' Suspected misspecification is represented with an induce distribution on the network
#' Currently only support cross-clusters contamination in cluster-randomized-trials (CRT)
#'
#' @param Z.obs Numeric vector of observed treatments
#' @param Y.obs Numeric vector of observed outcomes
#' @param theta_vec Vector of Theta parameters to run SA on
#' @param N_clusters Number of clusters present
#' @param N_each_cluster_vec Numeric vector with size N_clusters that indicate the number of units within each cluster
#' @param n.iter Number of SA iteration per Theta value
#' @param distance_matrix Optional between-clusters distance matrix; Deafult is NULL
#' @param distance.to.prob.func Function that convert distances and Theta value to probailities
#' @param Pz_function Function that generate treatment vectors (the experimental design function)
#' @param pz_func_args Parameters (as list) required by Pz_function
#' @param R Number of re-sampling from Pz_function to estimate the exposures probabilities
#'
#' @return Distribution of the causal effects as a function of theta (as data table)
#' @examples In the following example, the effects  tau(c11,c00), tau(c10,c00) are estimated
#' Z.obs <- Z_ber_clusters(N_clusters = 500, N_each_cluster_vec = rep(10,500), p = 0.5)
#' thres <- runif(5000)
#' base.Q <- generate_Q_matrix(N_clusters = 500,
#'                             within_prob_vec = rep(1,500),
#'                             between_prob_vec = 0.05)
#' a.mat <- as_adjacency_matrix(sample_sbm(n = 5000, pref.matrix = base.Q, block.sizes = rep(10,500), directed = FALSE))
#' expos <- generate_exposures_threshold(A = a.mat, Z = Z.obs, threshold = thres)
#' Y.obs <- generate_po(exposures = expos)
#' SensitivityAnalysis(Z.obs = Z.obs, Y.obs = Y.obs,
#'                     theta_vec = seq(0,0.1,0.05),
#'                     N_clusters = 500,
#'                     N_each_cluster_vec = rep(10,500),
#'                     n.iter = 100,
#'                     Pz_function = Z_ber_clusters,
#'                     pz_func_args = list(N_clusters = 500,
#'                                         N_each_cluster_vec = rep(10,500),
#'                                         p = 0.5),
#'                     R = 10^4)
#'
#'
#' @import igraph
#' @import data.table
#' @import gtools
#' @import parallel
#' @import mgcv
#'
SensitivityAnalysis <- function(Z.obs,
                                Y.obs,
                                theta_vec,
                                N_clusters,
                                N_each_cluster_vec,
                                n.iter,
                                distance_matrix = NULL,
                                distance.to.prob.func = exp_spatial_distance,
                                Pz_function,
                                pz_func_args,
                                R = 10^4){
  n <- length(Y.obs)
  # Run SA
  SA.results <- lapply(theta_vec, function(theta){
    data.table::rbindlist(lapply(seq(n.iter), function(m){

      # Generate block-Q matrix
      # Fixed Prob (=theta)
      if(is.null(distance_matrix)){
        Q.mat <- generate_Q_matrix(N_clusters = N_clusters,
                                   within_prob_vec = rep(1,N_clusters),
                                   between_prob_vec = theta)
      }
      # Prob by inverse of distance function (with parameter theta)
      else{
        dist_vec = distance_matrix[lower.tri(distance_matrix)]
        Q.mat <- generate_Q_matrix(N_clusters = N_clusters,
                                   within_prob_vec = rep(1, N_clusters),
                                   between_prob_vec = distance.to.prob.func(dist_vec, theta))
      }
      # Sample SBM
      adj.mat <- igraph::as_adjacency_matrix(igraph::sample_sbm(n = n,
                                                pref.matrix = Q.mat,
                                                block.sizes = N_each_cluster_vec,
                                                directed = FALSE))

      CE.estimate <- NMR_estimator(A.list = list(A=adj.mat),
                                   Z.obs = Z.obs,
                                   Y.obs = Y.obs,
                                   Pz_function = Pz_function,
                                   pz_func_args = pz_func_args,
                                   exposures_vec = c("c11","c00"),
                                   exposures_contrast = list(c("c11","c00")),
                                   exposure_func = generate_exposures_threshold,
                                   exposure_func_args = list(threshold = rep(0,n)),
                                   R = R)
      CE.estimate$iter <- m
      CE.estimate$theta <- theta
      CE.estimate
    }))
  })
  return(data.table::rbindlist(SA.results))
}

