
# source("R/Aux_functions.R")


#' Probabilistic bias analysis (PBA) for misspecification of the network interference structure
#'
#' @description
#' Compute the distribution of causal estimates obtained from perturbing the specified network with a suspected network deviation distribution and the relevant prior of choice.
#'
#' @details
#' Approximate Bayesian procedure that estimate causal estimates while accounting for systematic (deviations of the specified network) and random errors;
#' Suspected misspecification is represented by the induced distribution on the network;
#'
#' @param N_units Number of units (nodes) in the sample
#' @param N_iterations Number of PBA iterations (defaulted to 1000)
#' @param A.sp Baseline specified network
#' @param edges_prob_func A function that takes as input a parameter vector and specified network, and return a matrix with edge creation probabilities in the lower-triangle (see example)
#' @param prior_func_list Function of theta prior
#' @param prior_func_args_list list of parameters for the prior function
#' @param Z.obs Numeric vector of observed treatments
#' @param Y.obs Numeric vector of observed outcomes
#' @param X.obs Optional matrix or list of control variables that can be used in as inputs to ``edge_prob_func"
#' @param Pz_function Function that generate treatment vectors (the experimental design function)
#' @param pz_func_args Parameters (as list) required by Pz_function
#' @param exposures_vec Vector of exposures to compute
#' @param exposures_contrast List of causal contrast of interest
#' @param exposures_thresholds Optional thresholds for the exposure mapping
#' @param n.cores Optional number of cores for parallel computation (default is 1)
#'
#' @return N_iterations samples from the bias-corrected posterior of the point estimates (HT and Hajek) of causal effects with and without accounting for random error via normal approximation
#'
#' @examples
#' # In the following example, the effect  tau(c11,c00) is estimated with suspected deviation that randomly add missing edges with probability theta.
#' # The baseline network consists of 50 well separated clusters, but cross-cluster contamination is present.
#'
#' Z.obs <- Z_ber_clusters(N_clusters = 50, N_each_cluster_vec = rep(5,50), p = 0.5)
#' base.Q <- generate_Q_matrix(N_clusters = 50,
#'                             within_prob_vec = rep(1,250),
#'                             between_prob_vec = 0.05)
#' A.true <- as_adjacency_matrix(sample_sbm(n = 250, pref.matrix = base.Q, block.sizes = rep(5,50), directed = FALSE))
#' thres <- rep(0,250)
#' expos.true <- generate_exposures_threshold(A = A.true, Z = Z.obs, threshold = thres)
#' Y.obs <- generate_po(exposures = expos)
#'
# Q.sp <- generate_Q_matrix(N_clusters = 50,
#'                             within_prob_vec = rep(1,250),
#'                             between_prob_vec = 0)
#'A.sp <- as_adjacency_matrix(sample_sbm(n = 250, pref.matrix = Q.sp, block.sizes = rep(5,50), directed = FALSE))
#'
#'
#'clusters.id <- rep(seq(50),each=5)
#'
#' # Network deviation distribution assume to create edges between units in different clusters with prob. theta.
#'
#' contamination_edge_creation <- function(theta, A.sp, X){
#'  N <- dim(A.sp)[1]
#'  prob.mat <- matrix(NA, N, N)
#'  diag(prob.mat) <- 0
#'  for (i in seq(N-1)) {
#'    for (j in seq(i+1,N)) {
#'      # P(A_ij=1 | A^sp_ij = 0) = I(X[i]!=X[j])theta
#'      # P(A_ij=1 | A^sp_ij = 1) = 1
#'      prob.mat[j,i] <- ifelse(A.sp[j,i] == 1,
#'                              1,
#'                              (X[i] != X[j])*theta)
#'    }
#'  }
#'  return(prob.mat)
#'  }
#'
#'
#' # Run PBA with Uniform(0,0.005) prior for theta
#'
#' PBA.results <- PBA(N_units = 250,
#'                    N_iterations = 100,
#'                    A.sp = A.sp,
#'                    edges_prob_func = contamination_edge_creation,
#'                    prior_func_list = list(uniform),
#'                    prior_func_args_list = list(list(n=1,min=0,max=0.005)),
#'                    Z.obs = Z.obs,
#'                    Y.obs = Y.obs,
#'                    X.obs = clusters.id,
#'                    Pz_function = Z_ber_clusters,
#'                    pz_func_args = list(N_clusters = 50,
#'                                         N_each_cluster_vec = rep(5,50),
#'                                         p = 0.5),
#'                    exposures_vec = c("c11","c00"),
#'                    exposures_contrast = list(c("c11","c00"))
#'                    exposures_thresholds = thres)
#'
#' # Compare to baseline estimate
#'
#' baseline.estimate <- NMR_estimator(A.list = list(A.sp),
#'             Z.obs = Z.obs,
#'             Y.obs = Y.obs,
#'              Pz_function = Z_ber_clusters,
#'              pz_func_args = list(N_clusters = 50,
#'                                         N_each_cluster_vec = rep(5,50),
#'                                         p = 0.5),
#'             exposures_vec = c("c11","c00"),
#'            exposures_contrast = list(c("c11","c00")),
#'             exposure_func = generate_exposures_threshold,
#'             exposure_func_args = list(threshold = thres),
#'             R = 10^4)
#'
#' @import igraph
#' @import data.table
#' @import gtools
#' @import parallel
#' @import mgcv
#' @import Matrix
#' @export
#'
PBA <- function(N_units,
                N_iterations = 1000,
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
                exposures_thresholds = NULL,
                n.cores = 1){

  PBA.results <- mclapply(seq(N_iterations),function(m){
    CE.estimates <- Single_PBA_general_iteration(N_units = N_units,
                                                 A.sp = A.sp,
                                                 edges_prob_func = edges_prob_func,
                                                 prior_func_list = prior_func_list,
                                                 prior_func_args_list = prior_func_args_list,
                                                 Z.obs = Z.obs,
                                                 Y.obs = Y.obs,
                                                 X.obs = X.obs,
                                                 Pz_function = Pz_function,
                                                 pz_func_args = pz_func_args,
                                                 exposures_vec = exposures_vec,
                                                 exposures_contrast = exposures_contrast,
                                                 exposures_thresholds = exposures_thresholds)
    CE.estimates[,iter := m]
    CE.estimates
  },
  mc.cores = n.cores
  )
  return(rbindlist(PBA.results))
}
