library(ggplot2)
library(magrittr)

# >>> parameters >>>

## choose subject
#subject <- "K4B"
subject <- "K5B"
#subject <- "K6B"
#subject <- "K7B"
#subject <- "K8B"

#path_to_wes <- "C:/Users/kazik/outs/K4B/Strelka_S8934_FB_variants.vcf"
path_to_wes <- "C:/Users/kazik/outs/K5B/Strelka_S8934_2_FB_variants.vcf"
#path_to_wes <- "C:/Users/kazik/outs/K6B/Strelka_S13530_LN_vs_S13530_PBL_somatic_snvs.anno.vcf"
#path_to_wes <- "C:/Users/kazik/outs/K7B/Strelka_GS104753_LN_vs_GS104753_FB_somatic_snvs.anno.vcf"
#path_to_wes <- "C:/Users/kazik/outs/K8B/Strelka_GS104770_FB_variants.vcf"

save_path <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/germline/"
#save_path <- "../germline/'

#path_to_tumour <- "C:/Users/kazik/outs/K4B/Strelka_S8934_66747_FL_variants.anno.vcf"
path_to_tumour <- "C:/Users/kazik/outs/K5B/Strelka_S8934_1336260_FL_variants.anno.vcf"
#path_to_tumour <- "C:/Users/kazik/outs/K6B/Strelka_S13530_LN_vs_S13530_PBL_somatic_snvs.anno.vcf"
#path_to_tumour <- "C:/Users/kazik/outs/K7B/Strelka_GS104753_LN_vs_GS104753_FB_somatic_snvs.anno.vcf"
#path_to_tumour <- "C:/Users/kazik/outs/K8B/Strelka_GS104770_FB_variants.vcf"

# <<< parameters <<<


# >>> helper functions >>>
# <<< helper functions <<<


# >>> get normal positions >>>
wes_mutations <- read.delim(path_to_wes, 
                            header= FALSE, 
                            na.strings= ".", 
                            comment.char= "#")

#use only important germline mutations
wes_mutations <- wes_mutations[ wes_mutations$V7 == "PASS", c(1,2,4,5,10)]

#use only germline mutations of single nucleotides
ids_1 <- wes_mutations$V4 %>% sapply(function(n) is.element(n, c("A", "C", "G", "T")),
                                     USE.NAMES = FALSE)
ids_2 <- wes_mutations$V5 %>% sapply(function(n) is.element(n, c("A", "C", "G", "T")),
                                     USE.NAMES = FALSE)

ids <- ids_1 * ids_2

wes_mutations <- wes_mutations[ which(ids==1), ]
# <<< get normal positions <<<


# >>> get tumour positions >>>
tumour_reads <- read.delim(path_to_tumour, 
                           header= FALSE, 
                           na.strings= ".", 
                           comment.char= "#")

# use only important positions 
tumour_reads <- tumour_reads[ tumour_reads$V7 == "PASS", c(1,2,4,5,10)]

#use only germline mutations of single nucleotides
ids_1 <- tumour_reads$V4 %>% sapply(function(n) is.element(n, c("A", "C", "G", "T")),
                                     USE.NAMES = FALSE)
ids_2 <- tumour_reads$V5 %>% sapply(function(n) is.element(n, c("A", "C", "G", "T")),
                                     USE.NAMES = FALSE)

ids <- ids_1 * ids_2

tumour_reads <- tumour_reads[ which(ids==1), ]
# <<< get tumour positions <<<


# >>> join both selected positions >>>
all <- merge( wes_mutations, tumour_reads, by=c("V1", "V2") )

wes_mutations <- NULL
tumour_reads <- NULL

gc()
# <<< join both at selected positions <<<


# >>> get counts for tumour and normal >>>
# normal
nuc_counts <- all$V10.x %>% 
  strsplit( ':' ) %>% 
  sapply(function(entry) {entry[6]} )

alt_ref <- nuc_counts %>% 
  strsplit( ',' ) %>% 
  sapply(function(entry) {entry %>% as.integer})

all$AN <- alt_ref[1,]
all$BN <- alt_ref[2,]

#tumour
nuc_counts <- all$V10.y %>% 
  strsplit( ':' ) %>% 
  sapply(function(entry) {entry[6]} )

alt_ref <- nuc_counts %>% 
  strsplit( ',' ) %>% 
  sapply(function(entry) {entry %>% as.integer})

all$AT <- alt_ref[1,]
all$BT <- alt_ref[2,]
# <<< get counts for tumour and normal <<<

# >>> clean df >>>
all$V10.x <- NULL
all$V10.y <- NULL

colnames(all) <- c("CHR", "POS", 
                   "REF.N", "ALT.N", 
                   "REF.T", "ALT.T",
                   "AN", "BN",
                   "AT", "BT")
# <<< clean df <<<

# >>> save >>>
write.csv(all, file=paste0(save_path, "germline_counts_", subject, ".csv"))
# <<< save <<<
