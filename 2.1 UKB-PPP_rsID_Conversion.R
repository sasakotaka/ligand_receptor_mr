
## Load libraries

library(vroom)
library(tidyverse)
library(glue)
library(data.table)


## Load cytokine list

file <- 'Protein_Olink_Directory_Size_Chr28.txt'
protein_olink_category_list28 <- vroom(file, delim = "\t", col_name = TRUE)


## rsID conversion

for (j in 18){
  protein = protein_olink_category_list28$Protein[j]                                        ## protein = CXCL6 (j = 18)
  command1 <- paste0("file <- '/", protein, "/", protein, "_ppp_signif.txt'", sep = "")
  eval(parse(text = command1))
  snps_only_sumstat <- vroom(file, col_name = TRUE)
  snps_only_sumstat <- snps_only_sumstat %>%
    mutate(
      chrpos_id_1 = glue("{CHROM}:{GENPOS_h37}:{ALLELE0}:{ALLELE1}:imp:v1"),
      chrpos_id_2 = glue("{CHROM}:{GENPOS_h37}:{ALLELE1}:{ALLELE0}:imp:v1")
    )

  ## convert chr:pos to rsid
  rsid_chrpos_dir2 <- "SNP_RSID_maps"
  
  ## convert all the SNPs chr:pos to RSID, doing it per chromosomes
  rsid_snps_only_sumstat  <- data.frame()
  chr_to_run <- unique(snps_only_sumstat$CHROM)
  for (each_chr in chr_to_run){
    chr_info <- snps_only_sumstat %>% filter(CHROM == each_chr)
    chr_rsid_map <- read.table(
      file.path(
        rsid_chrpos_dir2, glue("rsid_map_chr{each_chr}.tsv")
      ), sep = "\t", header = TRUE
    ) 
    chr_rsid_map_list  <- setNames(
      chr_rsid_map$rsid, nm = chr_rsid_map$ID
    )
    chr_info$rsid <- NA
    chr_info$rsid <- ifelse(
      is.na(chr_info$rsid), chr_rsid_map_list[chr_info$ID], chr_info$rsid
    )
    rsid_snps_only_sumstat <- rbind(
      rsid_snps_only_sumstat,
      chr_info
    )
    command2 <- paste0("write.table(file='/", protein, "/", protein, "_ppp_signif_rsid.txt', rsid_snps_only_sumstat, sep = \"\t\", quote = FALSE, col.names = TRUE, row.names = FALSE)", sep = "")
    eval(parse(text = command2))
  }
}
