library(ggplot2)
library(magrittr)

# >>> parameters >>>

pvalue_cutoff <- .05
pval_str <- "005"

## choose subject
#subject <- "K4B"
#subject <- "K5B"
#subject <- "K6B"
subject <- "K7B"
#subject <- "K8B"

path_to_scrna <- "C:/Users/kazik/outs/K45678B.processed_somatic_snvs.ac"

#path_to_wes <- "C:/Users/kazik/outs/K4B/Strelka_S8934_66747_FL_vs_S8934_FB_somatic_snvs.anno.vcf"
#path_to_wes <- "C:/Users/kazik/outs/K5B/Strelka_S8934_1336260_FL_vs_S8934_FB_somatic_snvs.anno.vcf"
#path_to_wes <- "C:/Users/kazik/outs/K6B/Strelka_S13530_LN_vs_S13530_PBL_somatic_snvs.anno.vcf"
path_to_wes <- "C:/Users/kazik/outs/K7B/Strelka_GS104753_LN_vs_GS104753_FB_somatic_snvs.anno.vcf"
#path_to_wes <- "C:/Users/kazik/outs/K8B/Strelka_GS104770_LN_vs_GS104770_FB_somatic_snvs.anno.vcf"

save_path <- "C:/Users/kazik/FL_10X_2/350_CACTUS_extension/variants/"
#save_path <- "../variants/'

# <<< parameters <<<


# >>> helper functions >>>

pvalue_from_qual <- function(qual){10 ** (qual / (-10))}

# <<< helper functions <<<


# >>> read in scRNA >>>

all_scrna_counts <- read.delim( path_to_scrna )

# <<< read in scRNA <<<


# >>> significant wes mutations extraction >>>

wes_mutations <- read.delim(path_to_wes, 
                            header= FALSE, 
                            na.strings= ".", 
                            comment.char= "#")

qual_list <- wes_mutations$V8 %>% 
                strsplit( ';' ) %>% 
                sapply(function(entry) {entry[2]} )

qual_values <- qual_list %>% 
                strsplit( '=' ) %>% 
                sapply(function(entry) {entry[2] %>% as.integer} )

pvalues <- qual_values %>% sapply( pvalue_from_qual )

significant_wes_mutations <- wes_mutations[ pvalues < pvalue_cutoff , c(1,2,4,5)]

# <<< significant wes mutations extraction <<<


# >>> intersection of significant wes mutations and scrna >>>

scrna_counts <- all_scrna_counts[ all_scrna_counts$source == subject ,
                                  c("refPos", "refContig", "mixCnt", "altCnt")]
scrna_counts <- scrna_counts[ scrna_counts$mixCnt < 1 ,]
scrna_counts <- scrna_counts[ scrna_counts$altCnt > 0 , c("refPos", "refContig")]

scrna_by_chromosome <- split( scrna_counts, 
                              scrna_counts$refContig )
sig_wes_by_chromosome <- split( significant_wes_mutations,
                                significant_wes_mutations$V1 )

result <- data.frame("refContig"=NA, "refPos"=NA)

for ( chr in unique(significant_wes_mutations$V1) ){
  temp <- intersect(
    scrna_by_chromosome[[chr]]$refPos + 1,
    sig_wes_by_chromosome[[chr]]$V2
    )
  
  if( length(temp) ){
  result <- data.frame("refContig"=chr, "refPos"=temp) %>%
            rbind(result)
  }
}


num_wes <- dim( significant_wes_mutations )[1]
num_scrna <- unique( scrna_counts$refPos ) %>% length
num_intersection <- dim( result )[1]


saveRDS(result,
        file= stringr::str_c(save_path, subject, "/intersection_pvalue_", pval_str, ".rds")
        )


