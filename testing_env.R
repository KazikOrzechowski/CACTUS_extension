source('BCR_sim.R')
source('cactusx_clone_assignment.R')


for(pack in c('sp', 'fossil', 'ggplot2', 'doParallel')){
  if(!require(pack, lib.loc='libraries', character.only = TRUE)){
    install.packages(pack, lib='libraries', dependencies = TRUE)
    library(pack, lib.loc='libraries', character.only = TRUE)
  }
}


###############################################################################
#set hyperparameters
params <- read.csv('sim_params.csv')

#number of intermittent Gibbs sampling steps
i_G <- 3

max_iter <- 100000
buin_frac <- .5

new_sim <- TRUE
new_runs <- TRUE

n_runs <- 12
n_proc <- 12

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
    assignments <- cactus_clone_assignment(A = sim_A,
                                           D = sim_D,
                                           Omega = sim_Omega,
                                           BCR = sim_BCR,
                                           n_chain = n_runs,
                                           relax_rate_prior = c(0.5, 9.5),
                                           alpha_0=NULL,
                                           max_iter = max_iter,
                                           i_G= i_G,
                                           n_proc=n_proc
                                           )
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
    
    #pdf(file=paste0(save_path, str(run_), '_err_rate.pdf'), width=4, height=4)
    #print(ggplot() + geom_line(aes(x=1:max_iter, y=error_rate)))
    #dev.off()
    
    metrics$error_rates[run_] <- mean(error_rate[(buin_frac*max_iter) : max_iter])
    
    #calc and plot C matrix reconstruction. save mean and 0.1 quantile.
    C_prob <- dbinom(true_C, 1, assignments[[run_]]$C_prob)
    
    #pdf(file=paste0(save_path, str(run_), '_C_prob.pdf'), width=4, height=4)
    #print(ggplot() + geom_point(aes(x=1:(num_clones*num_variants), y=C_prob)))
    #dev.off()
    
    metrics$C_mean_aggre[run_] <- mean(C_prob)
    metrics$C_quant_aggre[run_] <- quantile(C_prob, probs=0.1)
    
    #calc and save assignment confidence
    x <- paste0('clone ', 1:num_clones)
    y <- 1:num_cells
    data <- expand.grid(X=x, Y=y)
    data$prob <- as.vector(t(assignments[[run_]]$prob_mat_j[order(assignments[[run_]]$t),]))
    
    #pdf(file=paste0(save_path, str(run_), '_ass_conf.pdf'), width=4, height=4)
    #print(ggplot(data) + geom_tile(aes(X, Y, fill=prob)))
    #dev.off()
    
    #calc and save B profile reconstruction error
    inferred_B <- assignments[[run_]]$B[assignments[[run_]]$t,,]
    
    av_abs_err <- sapply(1:num_cells, function(cell){mean(abs(inferred_B[cell,,]-true_B[cell,,]))})
    
    #pdf(file=paste0(save_path, str(run_), '_B_err.pdf'), width=4, height=4)
    #print(ggplot() + geom_point(aes(x=1:num_cells, y=av_abs_err)))
    #dev.off()
    
    metrics$B_profile_errs[run_] <- mean(av_abs_err)
    
    #plot and save data log-likelihood
    log_Lik <- assignments[[run_]]$logLik[2:max_iter]
    #log_Lik <- assignments[[run_]]$logLik[100*(1:(max_iter/100))]
    
    #pdf(file=paste0(save_path, str(run_), '_logLik.pdf'), width=4, height=4)
    #print(ggplot() + geom_line(aes(x=2:max_iter, y=log_Lik)))
    #ggplot() + geom_line(aes(x=1:(max_iter/100), y=log_Lik))
    #dev.off()
    
    metrics$log_likes[run_] <- mean(log_Lik[(buin_frac*max_iter) : max_iter])
    
    #calc adj rand index
    metrics$adj_rand_inds[run_] <- adj.rand.index(assignments[[run_]]$t, true_t)
  }
  
  saveRDS(metrics, file=paste0(save_path, '_metrics.RData'))
  
  #free memory for next simulation and runs
  assignments <- NULL
  gc()
} 

##################################################################################
