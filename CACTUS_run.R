source("C:/My_model/CACTUS/cactusx_helper_funcs.R")
source("C:/Users/kazik/Documents/RStudio/BCR_sim.R")

for(pack in c('fossil', 'ggplot2', 'doParallel')){
  if(!require(pack, character.only = TRUE)){
    install.packages(pack, lib='libraries', dependencies = TRUE)
    library(pack, lib.loc='libraries', character.only = TRUE)
  }
}

if(!require('igraph')){
  install.packages('igraph', lib='libraries')
  library('igraph', lib.loc='libraries')
}
if(!require('extraDistr')){
  install.packages('extraDistr', lib='libraries')
  library('extraDistr', lib.loc='libraries')
}
if(!require('matrixStats')){
  install.packages('matrixStats', lib='libraries')
  library('matrixStats', lib.loc='libraries')
}

if(!require('stableGR')){
  install.packages('stableGR', lib='libraries')
  library('stableGR', lib.loc='libraries')
}

loglike_RNA <- function(cells, clone, C, theta1, theta0, A, D){
  sum(dbinom(A[,cells],
             D[,cells], 
             prob=C[,clone]*theta1+(1-C[,clone])*theta0,
             log=TRUE))
}


cactus_clone_id_Gibbs <- function(A, D, Omega, BCR, Psi=NULL, 
                                  alpha_0=NULL, clusters_=NULL, t_=NULL,
                                  relax_C=TRUE, relax_rate_fixed=NULL, 
                                  relax_rate_prior=c(1, 19), keep_base_clone=TRUE,
                                  prior0=c(0.2, 99.8), prior1=c(0.45, 0.55),
                                  min_iter=5000, max_iter=75000, buin_frac=0.75,
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
    t_and_clust <- ini_clust(BCR, M, L)
    
    #t_ is cell to cluster assign
    t_ <- t_and_clust$t
    
    #clusters_ is the list of clusters containing their members
    clusters_ <- t_and_clust$clust
    
    Q <- length(clusters_)
    
    G = matrix(0, nrow=Q, ncol=M)
    
    for(cell in 1:M){
      G[,cell] <- sapply(1:Q, function(q){sum(abs(BCR[cell,,]-BCR[clusters_[[q]][1],,]))/2})
    }
    
    I <- extraDistr::rcat(Q, Psi)
    
    p_q <- matrix(0, nrow=Q, ncol=M)
    for(j in 1:M){
      temp  <- exp(-G[,j])
      p_q[,j] <- temp/sum(temp)
    }
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
  prob_all <- matrix(0, nrow = max_iter, ncol = M*K)
  
  theta0_all <- matrix(0, nrow = max_iter, ncol = 1)
  theta1_all <- matrix(0, nrow = max_iter, ncol = n_element)
  
  C_all <- matrix(0, nrow = max_iter, ncol = N*K)
  relax_rate_all <- matrix(0, nrow = max_iter, ncol = 1)
  
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
  
  ## Gibbs sampling
  for (it in 2:max_iter) {
    theta0 <- theta0_all[it - 1, 1]
    theta1 <- theta1_all[it - 1,  ]
    
    for(cell in 1:M){
      clust_probs <- sapply(1:Q, function(q){log(p_q[q, cell])+
                                             sum(dbinom(A[,cell],
                                                        D[,cell], 
                                                        prob=C[,I[q]]*theta1+(1-C[,I[q]])*theta0,
                                                        log=TRUE))})
      clust_probs <- clust_probs - max(clust_probs)
      t_[cell] <- extraDistr::rcat(1, exp(clust_probs))
    }
    for(clust in 1:Q){
      clusters_[[clust]] <- which(t_==clust, arr.ind=TRUE)
    }

    #perform non-conjugate split-merge moves
    prob_all[it, ] <- sample_I()
    
    assign_all_j[it, ] <- I[t_]
    
    #Update C and relax rate
    sample_c_and_rr()
    
    relax_rate_all[it] <- relax_rate
    C_all[it, ] <- C
    
    # Sample theta with assigned clones
    sample_theta()
    
    # Calculate logLikelihood
    if(TRUE){
      logLik_all[it] <- get_logLik(A1, B1, C, assign_all_j[it, ], 
                                   theta0_all[it,1], theta1_all[it,])
      #print(paste0(it/max_iter*100, '% of max iterations.'))
    }
    #Check convergence.
    if ((it >= min_iter) && (it %% 1000 == 0)) {
      # browser()
      Converged_all <- abs(Geweke_Z(prob_all[1:it, ])) <= 0.01
      # Converged_all <- abs(Geweke_Z(Config_all[1:it, ])) <= 2
      
      if (verbose) {
        cat(paste0(round(mean(Converged_all, na.rm = TRUE), 3) * 100, 
                   "% converged.\n"))
      }
      if (mean(Converged_all, na.rm = TRUE) > 0.995) {break}
    }
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
  
  prob_mat_j <- matrix(0, nrow = M, ncol = K)
  for(cell in 1:M){
    prob_mat_j[cell,] <- sapply(1:K, function(k){sum(assign_all_j[n_buin:it,cell]==k)}) / (it-n_buin+1)
  }
  
  logLik_post <- get_logLik(A1, B1, C_prob, prob_mat_j, theta0, theta1)
  DIC <- devianceIC(logLik_all[n_buin:it], logLik_post)
  
  
  return_list <- list("theta0" = theta0, "theta1" = theta1,
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
                      "prob_all"=prob_all)
  
  return_list
}


###############################################################################
#set hyperparameters
params <- read.csv('C:/My_model/CACTUS/sim_params.csv')

max_iter <- 100000
buin_frac <- .75

new_sim <- FALSE
new_runs <- TRUE

n_runs <- 2
n_proc <- 2
################################################################################
#start running simulations

for(r_num in 1:nrow(params)){
  this_sim <- params[r_num,]
  
  save_path <- this_sim$type
  av_reads <- this_sim$av_reads
  bcr_mut_prop <- this_sim$bcr_mut_prop
  mut_freq <- this_sim$mut_freq
  sim_alpha_0 <- this_sim$alpha_0
  
  if(this_sim$new_g){
    sim_g <- data.frame('A'=rpois(300,1/100)+1, 
                        'C'=rpois(300,10)+1,
                        'G'=rpois(300,1/100)+1,
                        'T'=rpois(300,1/100)+1)
  }else{
    sim_g <- data.frame('A'=rpois(300,1/100)+1, 
                        'C'=rpois(300,200)+1,
                        'G'=rpois(300,1/100)+1,
                        'T'=rpois(300,1/100)+1)
  }
  
  ################################################################################
  #data simulation
  if(new_sim){    
    simulated_data <- data_simulation(alpha_0=sim_alpha_0,
                                      mut_freq=mut_freq,
                                      av_reads=av_reads)
  }else{
    simulated_data <- readRDS(file=paste0(save_path, '_sim_data.RData'))
  }
  ###########################################################################
  #observed variables
  sim_A <- simulated_data$A
  sim_D <- simulated_data$D
  sim_Omega <- simulated_data$Omega
  sim_BCR <- simulated_data$BCR
  
  num_cells <- dim(sim_A)[2]
  num_clones <- dim(sim_Omega)[2]
  num_variants <- dim(sim_Omega)[1]
  
  ###########################################################################
  #hidden variables
  true_t <- simulated_data$t
  true_B <- simulated_data$B[true_t,,]
  true_I <- simulated_data$I
  true_C <- simulated_data$C
  true_t0 <- simulated_data$theta0
  true_t1 <- simulated_data$theta1
  
  ##############################################################################
  #run the model
  if(new_runs){
    registerDoParallel(cores=n_proc)
    
    n.eff <- get('n.eff')
    assignments <- foreach::foreach(i=1:n_runs, .export = ls(environment())) %dopar% {
                     cactus_clone_id_Gibbs(A = sim_A,
                                           D = sim_D,
                                           Omega = sim_Omega,
                                           BCR = sim_BCR,
                                           relax_rate_prior = c(0.5, 9.5),
                                           alpha_0=NULL,
                                           max_iter = max_iter)
    }
  }else{
    assignments <- readRDS(file=paste0(save_path, '_runs.RData'))
  }
  
  ###############################################################################
  #save results and plots
  saveRDS(simulated_data, file=paste0(save_path, '_sim_data.RData'))
  saveRDS(assignments, file=paste0(save_path, "_runs.RData"))
  
  metrics <- data.frame(sim_type=save_path,
                        run=1:n_runs,
                        error_rates=NA, 
                        C_mean_aggre=NA, 
                        C_quant_aggre=NA,
                        B_profile_errs=NA,
                        log_likes=NA,
                        adj_rand_inds=NA)
  
  for(run_ in 1:n_runs){
    #calc and plot error rate. save mean error
    error_rate <- sapply(1:max_iter, function(it){sum(true_I[true_t] != assignments[[run_]]$assign_all_j[it,]) / num_cells})
    
    pdf(file=paste0(save_path, str(run_), '_err_rate.pdf'), width=4, height=4)
    print(ggplot() + geom_line(aes(x=1:max_iter, y=error_rate)))
    dev.off()
    
    metrics$error_rates[run_] <- mean(error_rate[(buin_frac*max_iter) : max_iter])
    
    #calc and plot C matrix reconstruction. save mean and 0.1 quantile.
    C_prob <- dbinom(true_C, 1, assignments[[run_]]$C_prob)
    
    pdf(file=paste0(save_path, str(run_), '_C_prob.pdf'), width=4, height=4)
    print(ggplot() + geom_point(aes(x=1:(num_clones*num_variants), y=C_prob)))
    dev.off()
    
    metrics$C_mean_aggre[run_] <- mean(C_prob)
    metrics$C_quant_aggre[run_] <- quantile(C_prob, probs=0.1)
    
    #calc and save assignment confidence
    x <- paste0('clone ', 1:num_clones)
    y <- 1:num_cells
    data <- expand.grid(X=x, Y=y)
    data$prob <- as.vector(t(assignments[[run_]]$prob_mat_j[order(assignments[[run_]]$t),]))
    
    pdf(file=paste0(save_path, str(run_), '_ass_conf.pdf'), width=4, height=4)
    print(ggplot(data) + geom_tile(aes(X, Y, fill=prob)))
    dev.off()
    
    #plot and save data log-likelihood
    log_Lik <- assignments[[run_]]$logLik[2:max_iter]
    #log_Lik <- assignments[[run_]]$logLik[100*(1:(max_iter/100))]
    
    pdf(file=paste0(save_path, str(run_), '_logLik.pdf'), width=4, height=4)
    print(ggplot() + geom_line(aes(x=2:max_iter, y=log_Lik)))
    #ggplot() + geom_line(aes(x=1:(max_iter/100), y=log_Lik))
    dev.off()
    
    metrics$log_likes[run_] <- mean(log_Lik[(buin_frac*max_iter) : (max_iter-1)])
    
    #calc adj rand index
    metrics$adj_rand_inds[run_] <- adj.rand.index(assignments[[run_]]$t, true_t)
  }
  
  saveRDS(metrics, file=paste0(save_path, '_metrics.RData'))
  
  #free memory for next simulation and runs
  assignments <- NULL
  gc()
} 


