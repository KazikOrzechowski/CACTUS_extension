source("cactusx_helper_funcs.R")

cactusx_clone_assignment <- function(A, D, Config = NULL, BCR, n_clone = NULL, Psi = NULL, 
                     relax_Config = TRUE, relax_rate_fixed = NULL,
                     n_chain = 1, n_proc = 1, 
                     verbose = TRUE, relax_rate_prior = c(1,9),
                     alpha_0 = 0.1, ...) {

  ## check input data
  if (!(all(rownames(A) == rownames(D))))
    stop("Rownames for A and D are not identical.")
  if (!(all(colnames(A) == colnames(D))))
    stop("Colnames for A and D are not identical.")
  if (is.null(Config) && is.null(n_clone))
    stop("Config and n_clone can't be NULL together.")
 
  ## Match exome-seq and scRNA-seq data
  if (!any(rownames(D) %in% rownames(Config)))
    stop("No matches in variant names between Config and D arguments.")
  
  ## match variants
  common_vars <- intersect(rownames(Config), rownames(D))
  A <- A[common_vars,, drop = FALSE]
  D <- D[common_vars,, drop = FALSE]
  Config <- Config[common_vars,, drop = FALSE]
  if (verbose)
    message(length(common_vars), " variants used for cell assignment.")
  
  ## pass data to specific functions
  #doMC::registerDoMC(n_proc)
  `%dopar%` <- foreach::`%dopar%`
  
  #Do Gibbs sampling n_chain times independently
  ids_list <- foreach::foreach(ii = 1:n_chain) %dopar% {
    cactus_clone_id_Gibbs(A=A, D=D, Config=Config, BCR = BCR, 
                          Psi = Psi, alpha_0=alpha_0,
                          relax_Config = relax_Config, 
                          relax_rate_fixed = relax_rate_fixed,
                          relax_rate_prior = relax_rate_prior,
                          verbose = verbose)
  }
  
  #Average the results over n_chain models
  ids_out <- ids_list[[1]]
  ids_out$n_chain <- 1
  if (n_chain > 1) {
    for (ii in seq(2, n_chain)) {
      ids_out$n_chain <- ids_out$n_chain + 1
      idx <- colMatch(ids_out$prob, ids_list[[ii]]$prob, force = TRUE)
      ids_out$prob <- ids_out$prob + ids_list[[ii]]$prob[, idx]
      ids_out$relax_rate <- ids_out$relax_rate + ids_list[[ii]]$relax_rate
      ids_out$Config_prob <- (ids_out$Config_prob + 
                                ids_list[[ii]]$Config_prob[, idx])
    }
    ids_out$prob <- ids_out$prob / n_chain
    ids_out$relax_rate <- ids_out$relax_rate / n_chain
    ids_out$Config_prob <- ids_out$Config_prob / n_chain
  }
  return(ids_out)
}

#function describing the GIbbs sampling steps
cactus_clone_id_Gibbs <- function(A, D, Config, BCR, Psi=NULL, alpha_0=0.1,
                           relax_Config=TRUE, relax_rate_fixed=NULL, 
                           relax_rate_prior=c(1, 9), keep_base_clone=TRUE,
                           prior0=c(0.2, 99.8), prior1=c(0.45, 0.55),
                           min_iter=5000, max_iter=20000, buin_frac=0.5,
                           relabel=FALSE, verbose=TRUE) {

  if(is.null(Psi)){
    Psi <- rep(1/ncol(Config), ncol(Config))
    }
  if(dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] ||
     dim(A)[1] != dim(Config)[1] || dim(Config)[2] != length(Psi)){
      stop(paste0("A and D must have the same size;\n ",
                  "A and Config must have the same number of variants;\n",
                  "Config and Psi must have the same number of clones"))
  }
  
  ## preprocessing
  N <- dim(A)[1]             # number of variants
  M <- dim(A)[2]             # number of cells
  K <- dim(Config)[2]        # number of clones
  #Config is the matrix of true genotypes
  #Psi is the matrix of cluster->clone assignment
  ########################################################################################################
  BCR <- BCR[BCR$cell %in% colnames(A),] 
  
  if(nrow(BCR)!= M){
    stop(paste0("Input cells are different"))
  }
  
  Q <- length(unique(BCR$cluster))                      
  # number of clusters
  #clusters is the list of cluster representatives
  clusters <- unique(BCR$cluster)
  
  t_ <- list()
  for(q in 1:Q){
    t_[[q]] <- BCR$cell[BCR$cluster==clusters[q]]
  }
  #t is the list of cell assignments
  ########################################################################################################

  A[which(D == 0)] <- NA
  D[which(D == 0)] <- NA
  A[(D > 0) & is.na(A)] <- 0
  
  C1 <- Config
  C0 <- 1 - Config
  A1 <- A                  #number of alteration reads
  B1 <- D - A              #number of reference reads
  W_log <- sum(lchoose(D, A), na.rm = TRUE)  #log binomial coefficients
  
  A1[is.na(A1)] <- 0
  B1[is.na(B1)] <- 0
  
  #reads number list for each clone
  #these are lists with numbers of agreements/disagreements between A, D and C0, C1
  S1_list <- list()
  S2_list <- list()
  S3_list <- list()
  S4_list <- list()
  for (k in seq_len(K)) {
    S1_list[[k]] <- A1 * C0[, k]
    S2_list[[k]] <- B1 * C0[, k]
    S3_list[[k]] <- A1 * C1[, k]
    S4_list[[k]] <- B1 * C1[, k]
  }
  
  ## Prepare for sampling
  
  idx_vec <- seq_len(N)
  idx_mat <- seq_len(N*M)
  
  n_element <- length(idx_vec)
  
  if (is.null(dim(prior1)) && length(prior1) == 2) {
    #two variable to a matrix
    prior1 <- t(matrix(rep(prior1, n_element), nrow = 2))
  }
  if (!is.matrix(prior1)) {
    stop("prior1 need to be a matrix of n_element x 2")
  }
  
  prob_all   <- matrix(0, nrow = max_iter, ncol = Q*K)
  logLik_mat <- matrix(0, nrow = M, ncol = K)
  logLik_mat_q <- matrix(0, nrow = Q, ncol = K)
  logLik_all <- matrix(0, nrow = max_iter, ncol = 1)
  assign_all <- matrix(0, nrow = max_iter, ncol = Q)
  assign_all_j <- matrix(0, nrow = max_iter, ncol = M)
  theta0_all <- matrix(0, nrow = max_iter, ncol = 1)
  theta1_all <- matrix(0, nrow = max_iter, ncol = n_element)
  Config_all <- matrix(0, nrow = max_iter, ncol = N*K)
  relax_rate_all <- matrix(0, nrow = max_iter, ncol = 1)
  alpha_0_all <- matrix(0, nrow = max_iter, ncol = 1)
  
  relax_rate <- relax_rate_prior[1] / (relax_rate_prior[1] + 
                                             relax_rate_prior[2])
  
  Config_new <- Config
  Config_prior <- Config
  Config_prior[Config == 1] <- 1 - relax_rate
  Config_prior[Config == 0] <- relax_rate
  if (keep_base_clone) {
    Config_prior[, 1] <- Config[, 1]}
  Config_prior_oddlog <- log(Config_prior) - log(1 - Config_prior)
  Iden_mat <- matrix(0, nrow = M, ncol = K)

  
  ## Random initialization
  theta0_all[1,1] <- stats::rbeta(1, prior0[1], prior0[2])
  theta1_all[1, ] <- stats::rbeta(rep(1,n_element), prior1[,1], prior1[,2])
  theta1 <- matrix(NA, nrow = N, ncol = M)

  
  ## Set parent env of all called functions to this env
  environment(update_prob_mat) <- environment()
  environment(sample_cluster_to_clone_assign) <- environment()
  environment(sample_c_and_rr) <- environment()
  environment(sample_theta) <- environment()
  environment(sample_alpha_0) <- environment()
  
  ## Gibbs sampling
  for (it in 2:max_iter) {
    theta0 <- theta0_all[it - 1, 1]
    theta1[idx_mat] <- theta1_all[it - 1,  ]
    
    #Update prob_mat p_{j,q}
    update_prob_mat()
    prob_all[it, ] <- prob_mat
    
    # Sample assignment
    sample_cluster_to_clone_assign()
    
    #Update config and relax rate
    sample_c_and_rr()
    
    relax_rate_all[it] <- relax_rate
    Config_all[it, ] <- Config_new
    
    # Sample theta with assigned clones (do the same for gamma)
    sample_theta()

    # Sample Tj
    
    # Sample alpha_0
    sample_alpha_0()
    alpha_0_all[it, 1] <- alpha_0
    
    # Calculate logLikelihood
    logLik_all[it] <- get_logLik(A1, B1, Config_new, assign_all_j[it, ], 
                                 theta0_all[it,1], theta1_all[it,])
    
    #Check convergence.
    if ((it >= min_iter) && (it %% 100 == 0)) {
      Converged_all <- abs(Geweke_Z(prob_all[1:it, ])) <= 2
      if (verbose) {
        cat(paste0(round(mean(Converged_all, na.rm = TRUE), 3) * 100, 
                   "% converged.\n"))
      }
      if (mean(Converged_all, na.rm = TRUE) > 0.995) {break}
    }
  }
  print(paste("Converged in", it, "iterations."))
  
  ## Return values
  n_buin = ceiling(it * buin_frac)
  
  a <- A1[idx_mat]
  d <- A1[idx_mat] + B1[idx_mat]
  binom_pdf1 <- binom_pdf0 <- rep(0, n_element)
  for (i in seq(n_buin, it)) {
    binom_pdf1 <- binom_pdf1 + stats::dbinom(a, size = d,
                                             prob = theta1_all[i,])
    binom_pdf0 <- binom_pdf0 + stats::dbinom(a, size = d,
                                             prob = theta0_all[i])
  }
  prob_variant <- matrix(NA, nrow = N, ncol = M)
  prob_variant[idx_mat] <- binom_pdf1 / (binom_pdf1 + binom_pdf0)
  row.names(prob_variant) <- row.names(A)
  colnames(prob_variant) <- colnames(A)
  
  if (relabel) {
    col_idx_use <- seq(K)
    for (ii in seq(n_buin, it)) {
      mat1 <- matrix(prob_all[ii - 1, ], nrow = M)
      mat2 <- matrix(prob_all[ii, ], nrow = M)
      
      if (ncol(mat1) <= 5) {
        idx <- colMatch(mat1, mat2, force = TRUE) }
      else {
        idx <- colMatch(mat1, mat2, force = FALSE) }
      col_idx_use <- col_idx_use[idx]
      prob_all[ii, ] <- matrix(prob_all[ii, ], nrow = M)[, col_idx_use]
      Config_all[ii, ] <- matrix(Config_all[ii, ], nrow = N)[, col_idx_use]
    }
  }
  prob_mat <- matrix(colMeans(prob_all[n_buin:it, ]), nrow = Q)
  row.names(prob_mat) <- clusters #colnames(A)
  colnames(prob_mat) <- colnames(Config)
  
  Config_prob <- Config
  Config_prob[, ] <- colMeans(Config_all[n_buin:it, ])
  
  theta0 <- mean(theta0_all[n_buin:it, ])
  theta1[idx_mat] <- colMeans(as.matrix(theta1_all[n_buin:it, ]))
  
  alpha_0 <- mean(alpha_0_all[n_buin:it,])
  
  prob_mat_j <- matrix(0, nrow = M, ncol = K)
  for (q in seq_len(Q)) {
    prob_mat_j[match(t_[[q]], colnames(A1)),] <- prob_mat[q,] 
  }
  logLik_post <- get_logLik(A1, B1, Config_prob, prob_mat_j, theta0, theta1)
  DIC <- devianceIC(logLik_all[n_buin:it], logLik_post)
  
  
  return_list <- list("theta0" = theta0, "theta1" = theta1,
                      "alpha0" = alpha_0,
                      "theta0_all" = as.matrix(theta0_all[1:it, ]),
                      "theta1_all" = as.matrix(theta1_all[1:it, ]),
                      "element" = idx_mat, "logLik" = logLik_all[1:it],
                      "prob_all" = prob_all[1:it,],
                      "prob" = prob_mat, "prob_variant" = prob_variant,
                      "relax_rate" = mean(relax_rate_all[n_buin:it]),
                      "Config_prob" = Config_prob,
                      "Config_all" = Config_all[1:it, ],
                      "relax_rate_all" = relax_rate_all[1:it], "DIC"=DIC)
  return_list
}
