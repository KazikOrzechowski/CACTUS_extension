library(CODEX)
library(falconx)
library(magrittr)


#path_to_germline <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/germline/germline_counts_K4B.csv"
path_to_germline <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/germline/germline_counts_K5B.csv"
#path_to_germline <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/germline/germline_counts_K8B.csv"

#input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K4B/falconx_germline_input.csv"
input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K5B/falconx_germline_input.csv"
#input_save <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/falconx/K8B/falconx_germline_input.csv"


#####################################################################################
# Apply CODEX/CODEX2 to get total coverage bias
#####################################################################################
germline <- read.csv(path_to_germline)

reads <- germline[, c("AN", "BN", "AT", "BT")] %>% as.matrix

colnames(reads) <- c("AN", "BN", "AT", "BT")

# total read depth
Y <- matrix( nrow=nrow(reads),
             ncol=2 ) 

Y[, 1] <- reads[,1] + reads[,2]
Y[, 2] <- reads[,3] + reads[,4]


# For each chromosome
# Get GC content from a 50bp window centered at the SNP

pos <- as.numeric(germline[,'POS'])
ref <- IRanges(start=pos-25,end=pos+25)

gc <- rep(NA, nrow(reads))
for(chr in unique(germline$CHR)){
  ids <- germline$CHR == chr
  
  gc[ids] <- getgc(chr, ref[ids])  
}

# omit is.na(gc) in normalization

ids1 <- is.na(gc)
ids2 <- !ids1

# normalization
normObj <- normalize(Y[ids2,], gc[ids2], K=1)

Yhat <- Y
colnames( Yhat ) <- c( "sN", "sT" )
Yhat[ids2,] <- round( normObj$Yhat[[1]], 0 )

#
my_matrix <- germline[, 2:5] %>% as.matrix()
colnames( my_matrix ) <- c( "CHROM", "POS", "REF", "ALT" )

# save falconx input
cbind( my_matrix, reads, Yhat ) %>% write.csv(file=input_save,
                                              row.names = FALSE)
