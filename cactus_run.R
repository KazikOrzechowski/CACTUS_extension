source("cactusx_clone_assignment.R")

#set paths for data
treePath <- "~/RStudio/tree/S144/Z.rds"
wesPath <- "~/RStudio/WES/wes.rds"
scRNA_Path <- '~/RStudio/scRNA/ac.rds'
bcrPath <- '~/RStudio/scBCR/scBCR_GEX.rds'

Z <- readRDS(treePath)

xSample <- readRDS(wesPath)$S144

sc <- readRDS(scRNAPath)$S144[1:100,]
sample_sc <- sc[sc$MUTATION_ID %in% xSample$MUT_ID_minus, ]

bcr <- readRDS(bcrPath)$S144

preprocess <- FALSE
if(preprocess){
  A_clone <- data.frame(unique(sample_sc$cell),
                        row.names = unique(sample_sc$cell))
  A_clone <- cbind(A_clone,matrix(0, nrow(A_clone), length(xSample$MUT_ID_minus)))
  colnames(A_clone) <- c("Cell",xSample$MUT_ID_minus)
  
  D_clone <- A_clone
  
  for(read in 1:nrow(sample_sc)){
    d <- D_clone[D_clone$Cell==sample_sc[read,]$cell,][[sample_sc[read,]$MUTATION_ID]]
    D_clone[D_clone$Cell==sample_sc[read,]$cell,][[sample_sc[read,]$MUTATION_ID]] <- d+1
    if(sample_sc[read,]$refAllele!=sample_sc[read,]$base){
      a <- A_clone[A_clone$Cell==sample_sc[read,]$cell,][[sample_sc[read,]$MUTATION_ID]]
      A_clone[A_clone$Cell==sample_sc[read,]$cell,][[sample_sc[read,]$MUTATION_ID]] <- a+1
    }
  }
  saveRDS(A_clone,"~/RStudio/A_clone.rds")
  saveRDS(D_clone,"~/RStudio/D_clone.rds")
}else{
  A_clone <- readRDS("~/RStudio/A_clone.rds")
  D_clone <- readRDS("~/RStudio/D_clone.rds")
  rownames( A_clone) <- unique(sample_sc$cell)
  rownames( D_clone) <- unique(sample_sc$cell)
}

bcr_clusts <- bcr[,c(1,1)]
colnames(bcr_clusts) <- c('cell', 'cluster')

assignments <- cactus_clone_assignment(A = t(A_clone[,2:ncol(A_clone)]),
                                       D = t(D_clone[,2:ncol(D_clone)]),
                                       Config = Z[[1]],
                                       BCR = bcr_clusts,
                                       n_chain = 1,n_proc = 30,
                                       relax_rate_prior = c(0.5, 9.5))             
