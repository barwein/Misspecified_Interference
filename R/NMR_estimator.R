
# source("R/Aux_functions.R")


 #' Estimation of causal effects with the `network-misspecification-robust` estimator
 #'
 #' @description
 #' Return the estimated causal effects coupled with conservative variance using HT and Hajek NMR estimators
 #'
 #' @details
 #' Estimation of causal effects using a list of networks (represented as adjacency matrices);
 #' The estimator is unbiased if at least one of the networks is correctly specified
 #'
 #' @param A.list List of one or more adjacency matrices (in matrix format)
 #' @param Z.obs Numeric vector of observed treatments
 #' @param Y.obs Numeric vector of observed outcomes
 #' @param Pz_function Function that generate treatment vectors (the experimental design function)
 #' @param pz_func_args Parameters (as list) required by Pz_function
 #' @param exposure_vec Vector of characters of the relevant exposure values one wish to estimate the causal effects of
 #' @param exposures_contrast List of vectors, where each vector contains two exposure values (as characters)
 #' @param exposure_func Function that generate exposure values
 #' @param exposure_func_args Parameters (as list) required by exposure_func
 #' @param R Number of re-sampling from Pz_function to estimate the exposures probabilities
 #'
 #' @return Estimated causal effects (as data table)
 #' @examples
 #' a1 <- as_adjacency_matrix(sample_gnp(1000,0.08))
 #' a2 <- as_adjacency_matrix(sample_gnp(1000,0.05))
 #' Z.obs <- rbinom(1000,1,0.5)
 #' thres <- runif(1000)
 #' expos <- generate_exposures_threshold(A = a1, Z = Z.obs, threshold = thres)
 #' Y.obs <- generate_po(exposures = expos)
 #' exposure_vec <- c("c11","c00","c10","c01")
 #' exposures_contrast <- list(c("c11","c00"), c("c10","c00"))
 #'
 #' NMR_estimator(A.list = list(a1,a2),
 #'             Z.obs = Z.obs,
 #'             Y.obs = Y.obs,
 #'             Pz_function = Z_ber,
 #'             pz_func_args = list(n=1000,p=0.5),
 #'             exposures_vec = exposure_vec,
 #'             exposures_contrast = exposures_contrast,
 #'             exposure_func = generate_exposures_threshold,
 #'             exposure_func_args = list(threshold = thres),
 #'              R = 10^4)
 #'
 #'
 #' @import igraph
 #' @import data.table
 #' @import gtools
 #' @import parallel
 #' @import mgcv
 #'
NMR_estimator <- function(A.list,
                          Z.obs,
                          Y.obs,
                          Pz_function,
                          pz_func_args,
                          exposures_vec,
                          exposures_contrast,
                          exposure_func = generate_exposures_threshold,
                          exposure_func_args,
                          R = 10^4){

  Nv <- length(Y.obs)
  # Get prob. matrices of A.list networks
  prob.mat <- Get_prob_matrices_list(R = R,
                                     n = Nv,
                                     Pz_function = Pz_function,
                                     pz_func_args = pz_func_args,
                                     A.list = A.list,
                                     exposures_contrast = exposures_contrast,
                                     exposures_vec = exposures_vec,
                                     exposure_func = exposure_func,
                                     exposure_func_args = exposure_func_args)

  # Get observed exposures matrix
  observed.exposures <- exposures_under_network_list(Z = Z.obs,
                                                     A.list = A.list,
                                                     exposure_func = exposure_func,
                                                     param.list = exposure_func_args)

  # Estimate CE
  CE.dt <- rbindlist(MR_CE_estimator(Z.obs = Z.obs,
                                     Y.obs = Y.obs,
                                     A.list = A.list,
                                     exposures_contrast = exposures_contrast,
                                     exposures_vec = exposures_vec,
                                     Prob_matrices_list = prob.mat,
                                     exposure_func = exposure_func,
                                     exposure_func_args = exposure_func_args))

  CE.dt$ce_contrast <- sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])})
  return(CE.dt)
}
