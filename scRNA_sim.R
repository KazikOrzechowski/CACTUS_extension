source('BCR_sim.R')

#number of clones
K <- 3

#number of mutation spots
N <- 30

#relax_rate prior and simulation
relax_rate_prior <- c(1,9)

relax_rate <- rbeta(1, relax_rate_prior[1], relax_rate_prior[2])

#theta priors and simulation
prior0<-c(0.2, 99.8)
prior1<-c(0.45, 0.55)

theta0 <- rbeta(1, prior0[1], prior0[2])
theta1 <- rbeta(N, prior1[1], prior1[2])

#Psi is the prior distribution of cluster-clone assignment
Psi <- replicate(K, 1/K)

#Omega and C simulation
Omega <- array(rbinom(N*K, 1, .3), dim=c(N, K))
p <- abs(Omega-relax_rate)
C <- array(rbinom(N*K, 1, prob=p), dim=c(N,K))

#cluster-clone assignment
I <- 1:M
I[1:length(clusters_)] <- rcat(length(clusters_), Psi)

#All read counts and mutation reads simulation
D <- array(rpois(N*M, 1), dim=c(N, M))
p <- C[,I[t_]]*theta1 + (1-C[,I[t_]])*theta0
A <- array(rbinom(N*M, D, prob=p), dim=c(N,M))

