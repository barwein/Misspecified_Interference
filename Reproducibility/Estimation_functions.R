
library(igraph)
library(foreach)
library(doParallel)
library(data.table)
library(mgcv)
library(parallel)
library(gtools)

# install.packages("devtools")
# devtools::install_github("barwein/CIwithMIS")
library(CIwithMIS)

# Simulations -----------------------------------------------------------

# generate_exposures_threshold <- function(A, Z, threshold,deg=NULL){
#   # Compute num neighbors
#   if(is.null(deg)){
#     deg <- apply(A, 1, sum)
#   }
#   # Compute signs
#   indirect_sign <- as.numeric((A %*% Z) > threshold*deg)
#
#   # Save exposures
#   exposures <- vector(length = nrow(A))
#
#   exposures[Z*indirect_sign==1] <- "c11"
#   exposures[(1-Z)*indirect_sign==1] <- "c01"
#   exposures[Z*(1-indirect_sign)==1] <- "c10"
#   exposures[(1-Z)*(1-indirect_sign)==1] <- "c00"
#
#   return(exposures)
# }

generate_po <- function(exposures, base.noise = NULL){

  if(is.null(base.noise)){
    n <- length(exposures)
    y.00 <- runif(n = n, min = 0.5, max = 1.5)
  } else{
    y.00 <- base.noise
  }
  Y <- 2*y.00*(exposures=="c11") +
    1.25*y.00*(exposures=="c01") +
    1.5*y.00*(exposures=="c10") +
    1*y.00*(exposures=="c00")

  return(Y)
}


Z_ber <- function(n,p){
  rbinom(n = n, size = 1, prob = p)
}

# Z_ber_clusters <- function(N_clusters,
#                            N_each_cluster_vec,
#                            p){
#
#   treated_cluster_indicator = rep(0, N_clusters)
#   treated_cluster_idx = sample(x = seq(N_clusters),
#                                size = N_clusters %/% (1/p),
#                                replace = FALSE)
#   treated_cluster_indicator[treated_cluster_idx] <- 1
#   unit_level_treatment_vec = rep(treated_cluster_indicator, N_each_cluster_vec)
#   return(unit_level_treatment_vec)
# }

compute_exposure_prob <- function(Pz_function, pz_func_args, R, A, threshold){

  # Sample R Z-vector and save in a matrix (nv*R)
  Z_mat <- replicate(R, do.call(Pz_function, pz_func_args))

  # Convert treatments to exposures
  expos_mat <- apply(Z_mat, 2, function(x){
    generate_exposures_threshold(A = A,Z = x, threshold = threshold)
  })

  # Compute prob. using additive smoothing
  p_11 <- apply(expos_mat,MARGIN = 1, function(x){(sum(x=="c11")+1)/(R+1)})
  p_01 <- apply(expos_mat,MARGIN = 1, function(x){(sum(x=="c01")+1)/(R+1)})
  p_10 <- apply(expos_mat,MARGIN = 1, function(x){(sum(x=="c10")+1)/(R+1)})
  p_00 <- apply(expos_mat,MARGIN = 1, function(x){(sum(x=="c00")+1)/(R+1)})

  return(list("p11" = p_11, "p01" = p_01,
              "p10" = p_10, "p00" = p_00))

}

binom_prob <- function(a,b,p,n){
  return(pbinom(q = b,size = n, prob = p)-pbinom(q = a,size = n, prob = p))
}

compute_exposure_prob_analytical_ber <- function(A, p, threshold){

  Nv <- nrow(A)

  # Compute degrees vectors
  deg_vec <- apply(A, 2, sum)

  # Compute general threshold and round below
  general_threshold <- floor(deg_vec*threshold)

  # Compute exposure prob.

  p_11 <- apply(cbind(deg_vec, general_threshold), 1,
                function(x){p*binom_prob(a=x[2], b=x[1], p=p, n=x[1])})

  p_01 <- apply(cbind(deg_vec, general_threshold), 1,
                function(x){(1-p)*binom_prob(a=x[2], b=x[1], p=p, n=x[1])})

  p_10 <- apply(cbind(deg_vec, general_threshold), 1,
                function(x){p*binom_prob(a=-1, b=x[2], p=p, n=x[1])})

  p_00 <- apply(cbind(deg_vec, general_threshold), 1,
                function(x){(1-p)*binom_prob(a=-1, b=x[2], p=p, n=x[1])})

  # Return results
  return(list("p11" = p_11, "p01" = p_01,
              "p10" = p_10, "p00" = p_00))

}

compute_joint_exposure_prob <- function(Pz_function, nv, R, A1, A2, p, threshold){

  deg_A1 <- apply(A1, 1 , sum)
  deg_A2 <- apply(A2, 1 , sum)

  P_mat <- foreach(m = seq(R), .combine = "cbind") %dopar% {

    source("Simulations/Testing/Aux_functions.R")

    z <- Pz_function(n=nv,p=p)

    exposures_A1 <- generate_exposures_threshold(A = A1,
                                                 Z = z,
                                                 threshold = threshold,
                                                 deg = deg_A1)
    exposures_A2 <- generate_exposures_threshold(A = A2,
                                                 Z = z,
                                                 threshold = threshold,
                                                 deg = deg_A2)

    c11.e <- as.numeric(exposures_A1=="c11" & exposures_A2=="c11")
    c01.e <- as.numeric(exposures_A1=="c01" & exposures_A2=="c01")
    c10.e <- as.numeric(exposures_A1=="c10" & exposures_A2=="c10")
    c00.e <- as.numeric(exposures_A1=="c00" & exposures_A2=="c00")

    data.frame(c11.e, c01.e, c10.e, c00.e)

  }

  p_11 <- (rowSums(P_mat[,seq(1,4*R,4)]) + 1) / (R+1)
  p_01 <- (rowSums(P_mat[,seq(2,4*R+1,4)]) + 1) / (R+1)
  p_10 <- (rowSums(P_mat[,seq(3,4*R+2,4)]) + 1) / (R+1)
  p_00 <- (rowSums(P_mat[,seq(4,4*R+3,4)]) + 1) / (R+1)

  return(list("p11" = p_11, "p01" = p_01,
              "p10" = p_10, "p00" = p_00))
}


compute_true_exposure_prob <- function(network, p_z){

  deg_vec <- degree(network)

  true_p11 <- p_z*(1-((1-p_z)^deg_vec))
  true_p01 <- (1-p_z)*(1-((1-p_z)^deg_vec))
  true_p10 <- p_z*((1-p_z)^deg_vec)
  true_p00 <- (1-p_z)*((1-p_z)^deg_vec)

  return(list(true_p11 = true_p11,
              true_p01 = true_p01,
              true_p10 = true_p10,
              true_p00 = true_p00))
}

compute_exposure_misclassification <- function(true_exposure,
                                               noisy_exposure){
  m.c11 <- sum(noisy_exposure == "c11" & true_exposure != "c11")
  m.c01 <- sum(noisy_exposure == "c01" & true_exposure != "c01")
  m.c10 <- sum(noisy_exposure == "c10" & true_exposure != "c10")
  m.c00 <- sum(noisy_exposure == "c00" & true_exposure != "c00")
  return(c(m.c11 = m.c11,
           m.c01 = m.c01,
           m.c10 = m.c10,
           m.c00 = m.c00))
}

create_noisy_adj_mat <- function(true_adj_mat,
                                 n_v,
                                 alpha_,
                                 beta_){

  # Save only lower triangular matrix

  lower_tri_idx <- lower.tri(true_adj_mat)

  noisy_adj_mat <- true_adj_mat

  noisy_adj_mat[!lower_tri_idx] <- NA

  # Remove existing ties with prob. beta_
  # Add ties with prob. alpha_

  ind_ties <- which(noisy_adj_mat==1)
  ind_no_ties <- which(noisy_adj_mat==0)
  noisy_adj_mat[ind_ties] <- rbinom(length(ind_ties),1,1-beta_)
  noisy_adj_mat[ind_no_ties] <- rbinom(length(ind_no_ties),1,alpha_)

  # Compute the full noisy matrix from the lower tri one

  # diag is 0
  diag(noisy_adj_mat) <- 0

  # Force the upper tri be the same as the lower tri

  noisy_adj_mat <- as.matrix(Matrix::forceSymmetric(noisy_adj_mat,uplo="L"))

  return(noisy_adj_mat)
}


sample_network <- function(net_mod,
                           Nv,
                           p_v = NULL){
  if(net_mod == "PA"){ # Preferential attachment
    return(sample_pa(n = Nv, power = 1, directed = FALSE))
  }
  if(net_mod == "ER"){ # Erdős–Rényi random graph
    return(sample_gnp(n = Nv, p = p_v))
  }
  if(net_mod == "SW"){ # Small-world graph
    G <- sample_smallworld(dim = 1,
                           size = Nv,
                           nei = 2,
                           p = 0.1)
    return(simplify(G, remove.loops = TRUE))
  }
}

sample_multiple_networks <- function(net_models,
                                     Nv,
                                     p_v = NULL){
  network_list <- list()
  adj_mat_list <- list()
  for (n_net in seq(length(net_models))){
    cur_net <- sample_network(net_mod = net_models[n_net],
                              Nv = Nv,
                              p_v = p_v)
    cur_adj_mat <- as.matrix(as_adjacency_matrix(cur_net))

    network_list[[n_net]] <- cur_net
    adj_mat_list[[n_net]] <- cur_adj_mat
  }
  return(list(network_list = network_list,
              adj_mat_list = adj_mat_list))
}



n_joint_exposures <- function(expos_A, expos_obs){
  n.11 <- sum(expos_A=="c11" & expos_obs=="c11")
  n.01 <- sum(expos_A=="c01" & expos_obs=="c01")
  n.10 <- sum(expos_A=="c10" & expos_obs=="c10")
  n.00 <- sum(expos_A=="c00" & expos_obs=="c00")
  return(c(n.11, n.01, n.10, n.00))
}


mod_rbind <- function(dt1, dt2){
  rbindlist(list(dt1,dt2))
}

One_network_CE_estimation <- function(true_adj_mat,
                                      noisy_adj_mat,
                                      network_model,
                                      M,
                                      Nv,
                                      pz,
                                      R,
                                      param,
                                      # alpha_,
                                      # beta_,
                                      threshold,
                                      exposures_vec,
                                      exposures_contrast,
                                      base.po.noise = NULL,
                                      jaccard_idx){
  # Get prob. matrices for current censored network
  prob.mat <- Get_prob_matrices_list(R = R, n = Nv,
                                     Pz_function = Z_ber, pz_func_args = list(n=Nv,p=pz),
                                     A.list = list(noisy_adj_mat),
                                     exposures_contrast = exposures_contrast,
                                     exposures_vec = exposures_vec,
                                     exposure_func = generate_exposures_threshold,
                                     # threshold = threshold)
                                     exposure_func_args = list(threshold))
  # Estimate CE M times
  sim_results <- mclapply(seq(M), function(m){
    # Sample treatment vector
    Z <- Z_ber(n = Nv, p = pz)
    # Compute true exposures
    true_expos <- generate_exposures_threshold(A = true_adj_mat,
                                               Z = Z,
                                               threshold = threshold)
    # Compute observed outcomes (using TRUE adj mat)
    true_outcomes <- generate_po(true_expos, base.po.noise)
    # Compute exposures
    noisy_expos <- generate_exposures_threshold(A = noisy_adj_mat,
                                                Z = Z,
                                                threshold = threshold)

    CE_estimate <- rbindlist(MR_CE_estimator(Z.obs = Z,
                                             Y.obs = true_outcomes,
                                             A.list = list(noisy_adj_mat),
                                             exposures_contrast = exposures_contrast,
                                             exposures_vec = exposures_vec,
                                             Prob_matrices_list = prob.mat,
                                             exposure_func = generate_exposures_threshold,
                                             exposure_func_args = list(threshold),
                                             compute.variance = FALSE))[,1:2]

    # Compute number of exposure mis-classifications
    expos_misclass <- compute_exposure_misclassification(true_exposure = true_expos,
                                                         noisy_exposure = noisy_expos)

    # Update results DT
    CE_estimate[,`:=`(ce_contrast = sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])}),
                      network_model = network_model,
                      param = paste0(param, collapse = ", "),
                      iter = m,
                      n.expos.misclass = paste0(names(expos_misclass),"=",expos_misclass, collapse = "; "),
                      jac = jaccard_idx
    )]


    # Estimate exposure prob. (using NOISY adj mat)
    # expos_prob <- compute_exposure_prob(Pz_function = Z_ber,
    #                                     nv = Nv,
    #                                     R = R,
    #                                     A = noisy_adj_mat,
    #                                     p = pz,
    #                                     threshold = threshold)
    #
    # expos_prob <- compute_exposure_prob_analytical_ber(A = noisy_adj_mat,
    #                                                    p = pz,
    #                                                    threshold = threshold)
    # Estimate CE (using NOISY exposures)
    # cur_ce <- estimate_ce(expos_vec = noisy_expos,
    #                       P_expos = expos_prob,
    #                       Y = true_outcomes)

    # data.table(y_hat_ht = unlist(cur_ce[1:4]),
    #            y_hat_hajek = unlist(cur_ce[5:8]),
    #            exposure = c("c11","c01","c10","c00"),
    #            exposure_prop_true = c(table(true_expos)/Nv),
    #            exposure_prop_noisy = c(table(noisy_expos)/Nv),
    #            n_expos_misclass = unlist(expos_misclass),
    #            network_model = network_model,
    #            param = paste0(param, collapse = ", "),
    #            iter_ = m)
    CE_estimate
  },
  mc.cores = 1)
  return(rbindlist(sim_results))
}

Cluster_randomization_CE_estimation <- function(N_clusters,
                                                N_each_cluster_vec,
                                                theta_,
                                                Nv,
                                                pz_cluster,
                                                # Z,
                                                threshold,
                                                M,
                                                distance_matrix=NULL,
                                                base.po.noise = NULL){

  # Compute exposures assuming clusters network structure
  # cluster_expos <- generate_exposures_threshold(A = clusters_adj_mat,
  #                                               Z = Z)

  # Get SBM Q matrix
  if(is.null(distance_matrix)){
    cur_Q <- generate_Q_matrix(N_clusters = N_clusters,
                               within_prob_vec = rep(1, N_clusters),
                               between_prob_vec = rep(theta_, choose(N_clusters,2)))
  } else{
    # Get lower tri
    dist_vec = distance_matrix[lower.tri(distance_matrix)]
    cur_Q <- generate_Q_matrix(N_clusters = N_clusters,
                               within_prob_vec = rep(1, N_clusters),
                               between_prob_vec = exp(-dist_vec/(mean(dist_vec)*theta_)))
  }

  # Sample SBM
  cur_sbm <- sample_sbm(n = Nv,
                        pref.matrix = cur_Q,
                        block.sizes = N_each_cluster_vec)

  sbm_adj_mat <- as.matrix(as_adjacency_matrix(cur_sbm))

  # Replicate M time for current theta_ value
  sim_results <- lapply(seq(M), function(m){

    # Sample treatment vector (unit level vector but cluster randomization)
    cur_Z <- Z_ber_clusters(N_clusters = N_clusters,
                            N_each_cluster_vec = N_each_cluster_vec,
                            p = pz_cluster)
    # Compute exposures
    cur_expos <- generate_exposures_threshold(A = sbm_adj_mat,
                                              Z = cur_Z,
                                              threshold = threshold)

    # Compute observed outcomes using adj mat from SBM
    cur_outcomes <- generate_po(cur_expos, base.po.noise)


    # Estimate ATE using HT
    n = sum(N_each_cluster_vec)
    cur_ate = n^{-1}*2*(sum(cur_outcomes[cur_Z==1]) - sum(cur_outcomes[cur_Z==0]))


    expos_table = paste0(names(table(cur_expos)),"=",table(cur_expos), collapse = "; ")

    data.table(y_hat = cur_ate,
               y_true = 1, # Oracle ATE
               theta_ = theta_,
               expos_table = expos_table,
               m = m)
  })

  return(rbindlist(sim_results))
}


Bias_of_clusters_ATE <- function(Nv,
                                 N_clusters,
                                 N_each_cluster_vec,
                                 theta_seq,
                                 pz_cluster,
                                 threshold,
                                 M,
                                 distance_matrix=NULL,
                                 base.po.noise = NULL){


  full_results <- mclapply(theta_seq, function(theta_){

    # Estimate ATE for current theta_
    Cluster_randomization_CE_estimation(N_clusters = N_clusters,
                                        N_each_cluster_vec = N_each_cluster_vec,
                                        theta_ = theta_,
                                        Nv = Nv,
                                        pz_cluster = pz_cluster,
                                        threshold = threshold,
                                        M = M,
                                        distance_matrix = distance_matrix,
                                        base.po.noise = base.po.noise)

  },
  mc.cores = 1)
  # mc.cores = availableCores()/2)


  return(rbindlist(full_results))
}



create_noisy_adj_mat <- function(true_adj_mat,
                                 n_v,
                                 alpha_,
                                 beta_){

  # Save only lower triangular matrix

  lower_tri_idx <- lower.tri(true_adj_mat)

  noisy_adj_mat <- true_adj_mat

  noisy_adj_mat[!lower_tri_idx] <- NA

  # Remove existing ties with prob. beta_
  # Add ties with prob. alpha_

  ind_ties <- which(noisy_adj_mat==1)
  ind_no_ties <- which(noisy_adj_mat==0)
  noisy_adj_mat[ind_ties] <- rbinom(length(ind_ties),1,1-beta_)
  noisy_adj_mat[ind_no_ties] <- rbinom(length(ind_no_ties),1,alpha_)

  # Compute the full noisy matrix from the lower tri one

  # diag is 0
  diag(noisy_adj_mat) <- 0

  # Force the upper tri be the same as the lower tri

  noisy_adj_mat <- as.matrix(Matrix::forceSymmetric(noisy_adj_mat,uplo="L"))

  return(noisy_adj_mat)
}



create_censored_adj_mat <- function(true_adj_mat, K){
  censored.A <- true_adj_mat
  # Compute degrees
  deg.A <- rowSums(true_adj_mat)
  # Find which nodes with degree > K (and use permutation of it)
  to.censored.idx <- permute(which(deg.A > K))
  for (idx in to.censored.idx) {
    # Get current `idx` peers
    peers <- censored.A[,idx]
    n.peers <- sum(peers)
    if(n.peers > K){
      # sample n.peers-K units and remove the edge
      peers.ties.idx <- which(peers == 1)
      peers.to.remove <- sample(peers.ties.idx, n.peers - K, F)
      peers[peers.to.remove] <- 0
      # Update adj. matrix
      censored.A[,idx] <- peers
      censored.A[idx,] <- peers
    }
  }
  return(censored.A)
}



Bias_of_Noisy_Network_Simulation <- function(Nv,
                                             pz,
                                             p_v,
                                             net_model_name,
                                             # K_vec,
                                             beta_vec,
                                             threshold,
                                             adj_true_network,
                                             M,
                                             R,
                                             exposures_vec,
                                             exposures_contrast,
                                             base.po.noise = NULL){
  # Run the simulation
  full_results <- mclapply(beta_vec, function(beta){

    # Create noisy network out of censored adj mat
    if(beta>0){
      noisy_adj_mat <- create_noisy_adj_mat(true_adj_mat = adj_true_network,
                                            n_v =  Nv,
                                            alpha_ = beta/100,
                                            beta_ = beta)
    } else{
      noisy_adj_mat <- adj_true_network
    }

    jaccard_idx <- jaccard_edgeset_similarity(adj_true_network,noisy_adj_mat)

    One_network_CE_estimation(true_adj_mat = adj_true_network,
                              noisy_adj_mat = noisy_adj_mat,
                              network_model = net_model_name,
                              M = M,
                              Nv = Nv,
                              pz = pz,
                              R = R,
                              # param = K,
                              param = beta,
                              # alpha_ = alph_bet[1],
                              # beta_ = alph_bet[2],
                              threshold = threshold,
                              exposures_vec = exposures_vec,
                              exposures_contrast = exposures_contrast,
                              base.po.noise = base.po.noise,
                              jaccard_idx = jaccard_idx)


  },
  mc.cores = 1)
  # mc.cores = detectCores())


  return(rbindlist(full_results))
}

Bias_of_Censored_Network_Simulation <- function(Nv,
                                                pz,
                                                p_v,
                                                net_model_name,
                                                K_vec,
                                                # beta_vec,
                                                threshold,
                                                adj_true_network,
                                                M,
                                                R,
                                                exposures_vec,
                                                exposures_contrast,
                                                base.po.noise = NULL){
  # Run the simulation
  full_results <- mclapply(K_vec, function(K){

    # Create censored network
    if(K > 0){
      censored_adj_mat <- create_censored_adj_mat(true_adj_mat = adj_true_network,
                                                  K = K)
    } else{
      censored_adj_mat <- adj_true_network
    }

    jaccard_idx <- jaccard_edgeset_similarity(adj_true_network,censored_adj_mat)

    One_network_CE_estimation(true_adj_mat = adj_true_network,
                              noisy_adj_mat = censored_adj_mat,
                              network_model = net_model_name,
                              M = M,
                              Nv = Nv,
                              pz = pz,
                              R = R,
                              param = K,
                              # param = beta,
                              # alpha_ = alph_bet[1],
                              # beta_ = alph_bet[2],
                              threshold = threshold,
                              exposures_vec = exposures_vec,
                              exposures_contrast = exposures_contrast,
                              base.po.noise = base.po.noise,
                              jaccard_idx = jaccard_idx)


  },
  mc.cores = 1)
  # mc.cores = detectCores())


  return(rbindlist(full_results))
}



bias_variance_full_iter_in_setup <- function(R,
                                             k,
                                             iter,
                                             Pz_function,
                                             pz,
                                             true.adj.mat,
                                             A.list,
                                             adj.mat.names,
                                             with.adj,
                                             exposures_contrast, exposures_vec,
                                             threshold,
                                             base.po.noise=NULL){
  # For specific `k` value, A.list (subset of original), and other parameters,
  # Compute Prob. matrices and run `iter` iteration

  n <- nrow(A.list[[1]])

  # Need to compute probabilities matrices only once
  Prob_mat_list <- Get_prob_matrices_list(R = R, n = n,
                                          Pz_function = Pz_function,
                                          pz_func_args = list(n=n, p=pz),
                                          A.list = A.list,
                                          exposures_contrast = exposures_contrast,
                                          exposures_vec = exposures_vec,
                                          exposure_func = generate_exposures_threshold,
                                          # threshold = threshold)
                                          exposure_func_args = list(threshold))

  # `iter` times sample new Z, generate PO using true.adj.mat and estimate
  # using MR estimator and adj. matrices of A.list

  multiple_iter_results <- mclapply(seq(iter), function(m){
    # Generate data
    cur.Z <- Pz_function(n=n, p=pz)
    true.expos <- generate_exposures_threshold(A = true.adj.mat,
                                               Z = cur.Z,
                                               threshold = threshold)
    Y.obs <- generate_po(true.expos, base.po.noise)

    # Estimate CE using MR
    CE.dt <- rbindlist(MR_CE_estimator(Z.obs = cur.Z,
                                       Y.obs = Y.obs,
                                       A.list = A.list,
                                       exposures_contrast = exposures_contrast,
                                       exposures_vec = exposures_vec,
                                       Prob_matrices_list = Prob_mat_list,
                                       exposure_func = generate_exposures_threshold,
                                       exposure_func_args = list(threshold)))
    # Update results DT
    CE.dt[,`:=`(ce_contrast = sapply(exposures_contrast, function(x){paste0(x[1],"-",x[2])}),
                iter = m,
                K = k,
                with.true.adj = with.adj,
                adj.mat.used = paste0(adj.mat.names,collapse = ", ")
    )]

  },
  mc.cores = 1
  # mc.cores = availableCores()/2)
  # detectCores()/2
  )

  return(rbindlist(multiple_iter_results))

}

bias_variance_tradeoff_one_k_fullrun <- function(R,
                                                 k,
                                                 iter,
                                                 Pz_function, pz,
                                                 true.adj.mat,
                                                 true.adj.name,
                                                 A.list,
                                                 exposures_contrast, exposures_vec,
                                                 threshold,
                                                 base.po.noise = NULL){
  # Run bias-variance trade off simulation for all k-length adj. matrices combination

  names_combin <- combinations(n = length(A.list),
                               # r = k,
                               r = length(A.list),
                               v = names(A.list))

  one_k_results <- apply(names_combin, 1, function(cur.names){

    cur.a.lst <- A.list[cur.names]
    with.adj <- (true.adj.name %in% cur.names)

    bias_variance_full_iter_in_setup(R = R,
                                     k = k,
                                     iter = iter,
                                     Pz_function = Pz_function,
                                     pz = pz,
                                     true.adj.mat = true.adj.mat,
                                     A.list = cur.a.lst,
                                     adj.mat.names = cur.names,
                                     with.adj = with.adj,
                                     exposures_contrast = exposures_contrast,
                                     exposures_vec = exposures_vec,
                                     threshold = threshold,
                                     base.po.noise = base.po.noise)

  })

  return(rbindlist(one_k_results))
}

Bias_variance_tradeoff_in_MR <- function(R,
                                         K,
                                         iter,
                                         Pz_function, pz,
                                         true.adj.mat,
                                         true.adj.name,
                                         A.list,
                                         exposures_contrast, exposures_vec,
                                         threshold,
                                         base.po.noise = NULL){
  # Run bias-variance simulations for the MR estimator

  bias_variance_results <- lapply(seq(K), function(k){

    bias_variance_tradeoff_one_k_fullrun(R = R,
                                         k = k,
                                         iter = iter,
                                         Pz_function = Pz_function,
                                         pz = pz,
                                         true.adj.mat = true.adj.mat,
                                         true.adj.name = true.adj.name,
                                         A.list = A.list,
                                         exposures_contrast = exposures_contrast,
                                         exposures_vec = exposures_vec,
                                         threshold = threshold,
                                         base.po.noise = base.po.noise)
  })

  return(rbindlist(bias_variance_results))

}


jaccard_edgeset_similarity <- function(A1, A2) {
  # Transform to graph object
  G1 <- graph_from_adjacency_matrix(A1)
  G2 <- graph_from_adjacency_matrix(A2)
  # Compute intersection and union of edges
  inter <- length(E(G1 %s% G2))
  un <- length(E(G1 %u% G2))
  if (un == 0) {
    0
  } else {
    inter/un
  }
}


# Data analysis -----------------------------------------------------------

compute_prob_indicator_palluck <- function(R, n, Pz_function, pz_func_args, A.list, exposures_vec){
  # A function that returns a list of nXr indicator matrices (one for each exposures),
  # that indicate if unit i in sample r is exposed to some ck under all A \in A.list
  # The output is used for computation of exposures probabilities


  # Init a list of nXR indicator matrices
  Ind.mat.list <- vector("list",length(exposures_vec))
  for (i in seq(length(exposures_vec))) {
    assign(paste0("I.",exposures_vec[i]), matrix(NA, nrow = n, ncol = R))
    Ind.mat.list[[i]] <- get(paste0("I.",exposures_vec[i]))
  }
  names(Ind.mat.list) <- exposures_vec

  # Sample Z ~ PZ R times and update the indicators matrices
  for (r in seq(R)) {
    # cur.Z <- Pz_function(n=n, p=pz, schid=schid)
    cur.Z <- do.call(Pz_function, pz_func_args)
    cur.exposure.mat <- exposures_under_network_list_palluck(Z=cur.Z$Z.vec,
                                                             school.treated = cur.Z$school.treated,
                                                             A.list=A.list)
    for (ck in exposures_vec) {
      # Update indiactor matrices for each exposures and r
      Ind.mat.list[[ck]][,r] <- specific_exposure_under_network_list(exposures_mat=cur.exposure.mat,
                                                                     ck=ck)
    }
  }
  return(Ind.mat.list)
}


exposures_under_network_list_palluck <- function(Z, school.treated, A.list){
  # for m adj. matrices in A.list and n units, this function compute the nXm matrix of exposures under
  # each adj. matrix for treatment vector Z
  return(
    sapply(A.list, function(A){compute.exposures.palluck(Z, school.treated, A)})
  )}


treat.randomization.palluck <- function(n, pz=0.5, schid){
  trt.vector <- vector("numeric", length(schid))
  uniq.schools <- unique(schid)
  # get treated schools
  trt.schools <- sample(x = uniq.schools,size = round(length(uniq.schools)/2), replace = FALSE)
  # sample unit level treatment for each treated school
  for (schl in trt.schools) {
    trt.vector[schid == schl] <- rbinom(n=sum(schid==schl), size = 1, prob = pz)
  }
  is.treated.school <- schid %in% trt.schools
  return(list(Z.vec = trt.vector, school.treated = is.treated.school))
}


compute.exposures.palluck <- function(Z, school.treated, adj.mat){
  # Compute num neighbors
  # Compute signs
  indirect_sign <- as.numeric((adj.mat %*% Z) > 0)
  deg <- rowSums(adj.matrix.either)
  # Save exposures
  exposures <- vector(length = nrow(adj.mat))

  exposures[as.logical(Z*indirect_sign*school.treated)] <- "c111"
  exposures[as.logical(Z*(1-indirect_sign)*school.treated)] <- "c101"
  exposures[as.logical((1-Z)*indirect_sign*school.treated)] <- "c011"
  exposures[as.logical((1-Z)*(1-indirect_sign)*school.treated)] <- "c001"
  exposures[as.logical(1-school.treated)] <- "c000"
  exposures[deg==0] <- "ISOLATED"

  return(exposures)
}

MR_CE_estimator_palluck <- function(Z.obs,school.treated ,Y.obs, A.list,
                                    exposures_contrast, exposures_vec,
                                    Prob_matrices_list,
                                    threshold,
                                    estimate.n = FALSE){

  expos.obs.mat <- exposures_under_network_list_palluck(Z = Z.obs,school.treated = school.treated,
                                                        A.list = A.list)

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
                                 estimate.n = estimate.n)
    ce.results[[cntrast_name]] <- ce_estimate
  }
  return(ce.results)
}




