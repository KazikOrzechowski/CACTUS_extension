source("cactusx_helper_funcs.R")
source('sampling_t.R')

for(pack in c('igraph', 'extraDistr', 'matrixStats')){
  if(!require(pack, lib.loc='libraries', character.only=TRUE)){
    install.packages(pack, lib='libraries')
    library(pack, lib.loc='libraries', character.only=TRUE)
  }
}

cactus_clone_assignment <- function(A, D, Omega = NULL, BCR, 
                                    n_clone = NULL, Psi = NULL, 
                                    relax_C = TRUE, relax_rate_fixed = NULL,
                                    n_chain = 1,  
                                    verbose = TRUE, relax_rate_prior = c(1,19),
                                    alpha_0 = NULL, max_iter=20000,
                                    i_G = 3, ...){

  ## check input data
  if (!(all(rownames(A) == rownames(D))))
    stop("Rownames for A and D are not identical.")
  if (!(all(colnames(A) == colnames(D))))
    stop("Colnames for A and D are not identical.")
  if (is.null(Omega) && is.null(n_clone))
    stop("C and n_clone can't be NULL together.")
 
  ## Match exome-seq and scRNA-seq data
  #if (!any(rownames(D) %in% rownames(Omega)))
  #  stop("No matches in variant names between C and D arguments.")
  
  ## match variants
  #common_vars <- intersect(rownames(Omega), rownames(D))
  #A <- A[common_vars,, drop = FALSE]
  #D <- D[common_vars,, drop = FALSE]
  #Omega <- Omega[common_vars,, drop = FALSE]
  #if (verbose)
  #  message(length(common_vars), " variants used for cell assignment.")
  
  ## pass data to specific functions
  registerDoParallel(cores=6)
  
  cactus_clone_id_Gibbs <- get('cactus_clone_id_Gibbs')
  #Do Gibbs sampling n_chain times independently
  ids_list <- foreach::foreach(i=1:n_chain, .export = ls(environment())) %dopar% {
    source("cactusx_helper_funcs.R")
    source('sampling_t.R')
    cactus_clone_id_Gibbs(A=A, D=D, Omega=Omega, BCR=BCR, 
                          Psi = Psi, alpha_0=alpha_0,
                          relax_C = relax_C, 
                          relax_rate_fixed = relax_rate_fixed,
                          relax_rate_prior = relax_rate_prior,
                          max_iter = max_iter,
                          i_G = i_G,
                          verbose = verbose)
  }
  #Average the results over n_chain models
  # ids_out <- ids_list[[1]]
  # ids_out$n_chain <- 1
  # if (n_chain > 1) {
  #   for (ii in seq(2, n_chain)) {
  #     ids_out$n_chain <- ids_out$n_chain + 1
  #     ids_out$relax_rate <- ids_out$relax_rate + ids_list[[ii]]$relax_rate
  #     ids_out$C_prob <- (ids_out$C_prob + 
  #                               ids_list[[ii]]$C_prob)
  #   }
  #   ids_out$relax_rate <- ids_out$relax_rate / n_chain
  #   ids_out$C_prob <- ids_out$C_prob / n_chain
  # }
  
  ids_list
}


#function describing the Gibbs sampling steps
cactus_clone_id_Gibbs <- function(A, D, Omega, BCR, Psi=NULL, 
                                  alpha_0=NULL, clusters_=NULL, t_=NULL,
                                  relax_C=TRUE, relax_rate_fixed=NULL, 
                                  relax_rate_prior=c(1, 19), keep_base_clone=TRUE,
                                  prior0=c(0.2, 99.8), prior1=c(0.45, 0.55),
                                  min_iter=5000, max_iter=40000, buin_frac=0.5,
                                  i_G = i_G,
                                  relabel=FALSE, verbose=TRUE) {
  
  if(is.null(Psi)){
    Psi <- rep(1/ncol(Omega), ncol(Omega))
  }
  if(dim(A)[1] != dim(D)[1] || dim(A)[2] != dim(D)[2] ||
     dim(A)[1] != dim(Omega)[1] || dim(Omega)[2] != length(Psi)){
    stop(paste0("A and D must have the same size;\n ",
                "A and Omega must have the same number of variants;\n",
                "Omega and Psi must have the same number of clones"))
  }
  
  ## preprocessing
  N <- dim(A)[1]             # number of variants
  M <- dim(A)[2]             # number of cells
  K <- dim(Omega)[2]         # number of clones
  L <- dim(BCR)[2]           # length of BCR seq
  #Omega is the matrix of inferred mutation profiles
  #Psi is prior of cluster->clone assignment
  ########################################################################################################
  if(is.null(clusters_) || is.null(t_)){
    #print('Initial clustering.')
    
    g <- matrix(1, nrow=L, ncol=4)
    
    t_and_clust <- ini_clust(BCR, M, L)
    
    #t_ is cell to cluster assign
    t_ <- t_and_clust$t
    
    #clusters_ is the list of clusters containing their members
    clusters_ <- t_and_clust$clust
    
    I <- extraDistr::rcat(length(clusters_), Psi)
    
    B <- array(0, dim=c(M, L, 4))
    
    for(i_ in 1:length(clusters_)){
      cl <- clusters_[[i_]]
      
      if(length(cl)>1){
        up_prior <- as.matrix(g) + apply(BCR[cl,,], c(2,3), sum)
      }
      if(length(cl)==1){
        up_prior <- as.matrix(g) + BCR[cl,,]
      }
      
      B[i_,,] <- t(sapply(1:L, function(l){igraph::sample_dirichlet(1, up_prior[l,])}))
    }
    
    alpha_0 <- 100
  }
  ########################################################################################################
  
  #A[which(D == 0)] <- NA
  #D[which(D == 0)] <- NA
  #A[(D > 0) & is.na(A)] <- 0
  
  C1 <- Omega
  C0 <- 1 - Omega
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
  
  logLik_all <- matrix(0, nrow = max_iter, ncol = 1)
  assign_all_j <- matrix(0, nrow = max_iter, ncol = M)
  
  theta0_all <- matrix(0, nrow = max_iter, ncol = 1)
  theta1_all <- matrix(0, nrow = max_iter, ncol = n_element)
  
  C_all <- matrix(0, nrow = max_iter, ncol = N*K)
  relax_rate_all <- matrix(0, nrow = max_iter, ncol = 1)
  
  alpha_0_all <- matrix(0, nrow = max_iter, ncol = 1)
  
  relax_rate <- relax_rate_prior[1] / (relax_rate_prior[1] + 
                                         relax_rate_prior[2])
  
  C <- Omega
  C_prior <- Omega
  C_prior[Omega == 1] <- 1 - relax_rate
  C_prior[Omega == 0] <- relax_rate
  C_prior_oddlog <- log(C_prior) - log(1 - C_prior)
  Iden_mat <- matrix(0, nrow = M, ncol = K)
  
  
  ## Random initialization
  theta0_all[1,1] <- stats::rbeta(1, prior0[1], prior0[2])
  theta1_all[1, ] <- stats::rbeta(rep(1,n_element), prior1[,1], prior1[,2])
  
  
  ## Set parent env of all called functions to this env
  environment(sample_I) <- environment()
  environment(sample_c_and_rr) <- environment()
  environment(sample_theta) <- environment()
  environment(sample_alpha_0) <- environment()
  environment(try_merge) <- environment()
  environment(try_split) <- environment()
  
  ## Gibbs sampling
  for (it in 2:max_iter) {
    theta0 <- theta0_all[it - 1, 1]
    theta1 <- theta1_all[it - 1,  ]
    
    #perform non-conjugate split-merge moves
    sample_I()
    try_merge()
    try_split()
    
    assign_all_j[it, ] <- I[t_]
    
    #Update C and relax rate
    sample_c_and_rr()
    
    relax_rate_all[it] <- relax_rate
    C_all[it, ] <- C
    
    # Sample theta with assigned clones
    sample_theta()
    
    # Sample alpha_0
    sample_alpha_0(alpha_prior=c(10, 0.1))
    
    alpha_0_all[it, 1] <- alpha_0
    
    # Calculate logLikelihood
    if(TRUE){
      logLik_all[it] <- get_logLik(A1, B1, C, assign_all_j[it, ], 
                                   theta0_all[it,1], theta1_all[it,]) +
        get_logLik_BCR(BCR, B, clusters_, t_)
      #print(paste0(it/max_iter*100, '% of max iterations.'))
    }
    #Check convergence.
    #if ((it >= min_iter) && (it %% 100 == 0)) {
    #  Converged_all <- abs(Geweke_Z(prob_all[1:it, ])) <= 2
    #  if (verbose) {
    #    cat(paste0(round(mean(Converged_all, na.rm = TRUE), 3) * 100, 
    #               "% converged.\n"))
    #  }
    #  if (mean(Converged_all, na.rm = TRUE) > 0.995) {break}
    #}
  }
  #print(paste("Converged in", it, "iterations."))
  
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
  
  C_prob <- matrix(colMeans(C_all[n_buin:it, ]), nrow=N, ncol=K)
  
  theta0 <- mean(theta0_all[n_buin:it, ])
  theta1[idx_mat] <- colMeans(as.matrix(theta1_all[n_buin:it, ]))
  
  alpha_0 <- mean(alpha_0_all[n_buin:it,])
  
  prob_mat_j <- matrix(0, nrow = M, ncol = K)
  for(cell in 1:M){
    prob_mat_j[cell,] <- sapply(1:K, function(k){sum(assign_all_j[n_buin:it,cell]==k)}) / (it-n_buin+1)
  }
  
  logLik_post <- get_logLik(A1, B1, C_prob, prob_mat_j, theta0, theta1)
  DIC <- devianceIC(logLik_all[n_buin:it], logLik_post)
  
  
  return_list <- list("theta0" = theta0, "theta1" = theta1,
                      "alpha0" = alpha_0_all,
                      "theta0_all" = as.matrix(theta0_all[1:it, ]),
                      "theta1_all" = as.matrix(theta1_all[1:it, ]),
                      "element" = idx_mat, "logLik" = logLik_all[1:it],
                      "prob_variant" = prob_variant,
                      "assign_all_j" = assign_all_j,
                      "prob_mat_j"= prob_mat_j,
                      "relax_rate" = mean(relax_rate_all[n_buin:it]),
                      "C_prob" = C_prob,
                      "C_all" = C_all[1:it, ],
                      "relax_rate_all" = relax_rate_all[1:it], "DIC"=DIC,
                      "clusters" = clusters_,
                      "t"= t_,
                      "I"= I,
                      "B"=B)
  
  return_list
}



