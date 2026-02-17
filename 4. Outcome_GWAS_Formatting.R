
## Load libraries

library('vroom')
library('glue')
library('TwoSampleMR')


## Load harmonized GWAS data from GWAS Catalog

#CRP
outco <- vroom('GCST90179146.h.txt')
outcome <- "CRP"

#RA
outco <- vroom('GCST90476230.h.txt')
outco$beta <- log(outco$odds_ratio)
outco$standard_error <- (log(outco$ci_upper) - log(outco$ci_lower)) / (2 * 1.96)
outcome <- "RA"

#Neurtophil count
outco <- vroom('32888493-GCST90002351-EFO_0004833.h.txt')
outcome <- "neutrophil"

#Monocyte count
outco <- vroom('32888493-GCST90002340-EFO_0005091.h.txt')
outcome <- "monocyte"

#Lymphocyte count
outco <- vroom('32888493-GCST90002316-EFO_0004587.h.txt')
outcome <- "lymphocyte"


## Format the harmonized GWAS data from GWAS Catalog

outco <- outco %>%
      mutate(
        chr_genpos = glue("{chromosome}_{base_pair_location}"),
      )    
outco_dat <- format_data(outco, header = TRUE, type = "outcome", 
                         beta_col = "beta",
                         se_col = "standard_error",
                         effect_allele_col = "effect_allele",
                         eaf_col = "effect_allele_frequency",
                         other_allele_col = "other_allele",
                         pval_col = "p_value",
                         chr_col = "chromosome",
                         pos_col = "base_pair_location",
                         snp_col = "chr_genpos",
                         log_pval = FALSE)
command <- paste0("write.table('outco_dat_", outcome, "_formatted.txt', outco_dat, sep = \"\t\", quote = FALSE, col.names = TRUE, row.names = FALSE)", sep = "")
eval(parse(text = command))


## Load GWAS data from FinnGen

#Bladder cancer
outco <- vroom('finngen_R12_C3_BLADDER_EXALLC_UKBB_meta.txt')
outcome <- "BladderCancer"


## Format the GWAS data from FinnGen

outco <- outco %>%
  mutate(
    chr_genpos = glue("{CHR}_{POS}"),
  )
outco_dat <- format_data(outco, header = TRUE, type = "outcome", 
                         beta_col = "all_inv_var_meta_beta", 
                         se_col = "all_inv_var_meta_sebeta", 
                         effect_allele_col = "ALT", 
                         other_allele_col = "REF", 
                         pval_col = "all_inv_var_meta_p", 
                         chr_col = "CHR", 
                         pos_col = "POS", 
                         snp_col = "chr_genpos", 
                         log_pval = FALSE)
command <- paste0("write.table('outco_dat_", outcome, "_formatted.txt', outco_dat, sep = \"\t\", quote = FALSE, col.names = TRUE, row.names = FALSE)", sep = "")
eval(parse(text = command))
