
## Load libraries

library('vroom')
library('dplyr')
library('tidyverse')
library('glue')
library('TwoSampleMR')


## Load data

file <- 'Protein_Olink_Directory_Size_Chr28.txt'
protein_list28 <- vroom(file, col_name = TRUE)

file <- 'outco_dat_CRP_formatted.txt'
outco_CRP <- vroom(file, col_names = TRUE)

file <- 'outco_dat_RA_formatted.txt'
outco_RA <- vroom(file, col_names = TRUE)

wbc_list <- c("CRP", "RA2024")
wbc_N <- c(575531, effective_n(10479, 436161))


## Perform MR

j <- 29
for (i in 1:2) {
  wbc = wbc_list[i]
  protein = protein_list28$Protein[j]
  protein_chr = protein_list28$Chr[j]
  protein_tss = protein_list28$TSS[j]
  
  command1 <- paste0("outco <- outco_", wbc, "
                    file <- '", protein, "_ppp_signif_extract.txt'", sep = "")
  eval(parse(text = command1))
  
  ligand_pqtl <- vroom(file, col_names = TRUE)
  ligand_pqtl <- ligand_pqtl %>%
    mutate(
      chr_genpos = glue("{CHROM}_{GENPOS}"),
    )
  ligand_pqtl <- ligand_pqtl %>%
    mutate(MHC = case_when(
      ligand_pqtl$CHROM == 6 & ligand_pqtl$GENPOS > 28010120 & ligand_pqtl$GENPOS < 33980577 ~ "YES",
      TRUE ~ "NO"))
  ligand_pqtl <- ligand_pqtl %>%
    mutate(CIS_TRANS = case_when(
      ligand_pqtl$CHROM == protein_chr & ligand_pqtl$GENPOS > protein_tss - 1000000 & ligand_pqtl$GENPOS < protein_tss + 1000000 ~ "cis",
      TRUE ~ "trans"))
  all <- filter(ligand_pqtl, MHC == "NO")
  trans <- filter(all, CIS_TRANS == "trans")
  receptor <- filter(trans, REC == "rec")
  non_receptor <- filter(trans, REC == "FALSE")
  
  
  ## MR using all pQTLs
  expo_dat_all <- format_data(all, type = "exposure", snp_col = "chr_genpos", beta_col = "BETA", se_col = "SE",
                              eaf_col = "EAF", effect_allele_col = "effect_allele", other_allele_col = "other_allele", 
                              pval_col = "log10P", chr_col = "chr", pos_col = "pos", log_pval = TRUE)
  dat_all <- harmonise_data(exposure_dat = expo_dat_all, outcome_dat = outco, action = 1)
  dat_all$samplesize.exposure <- protein_list28$N[j]
  dat_all$samplesize.outcome <- wbc_N[i]
  dat_all <- steiger_filtering(dat_all)
  dat_all$Rsq <- (dat_all$beta.exposure / dat_all$se.exposure) ^2 / ((dat_all$beta.exposure / dat_all$se.exposure) ^2 + dat_all$samplesize.exposure - 2)
  dat_all$F_stat <- dat_all$Rsq * (dat_all$samplesize.exposure - 2)/(1 - dat_all$Rsq)
  command2 <- paste0("write.table('mr_data_", wbc, "_", protein, ".txt', dat_all, sep = \"\t\", quote = FALSE, col.names = TRUE, row.names = FALSE)", sep = "")
  eval(parse(text = command2))
  
  res_mr_all <- mr(dat_all, method_list=c("mr_ivw", "mr_wald_ratio", "mr_egger_regression", "mr_weighted_median"))
  res_mr_all_odds <- generate_odds_ratios(res_mr_all)
  if (nrow(filter(dat_all, remove == "FALSE"))  > 2) {
    het <- mr_heterogeneity(dat_all, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  }
  if (nrow(filter(dat_all, remove == "FALSE"))  == 2) {
    het <- mr_heterogeneity(dat_all, method_list = c("mr_ivw", "mr_egger_regression"))
  }
  het$Isq <- (het$Q - het$Q_df) * 100 / het$Q
  res_mr_all_odds_het <- left_join(res_mr_all_odds, het, by = 'method')
  plt <- mr_pleiotropy_test(dat_all) 
  res_mr_all_odds_het_plt <- left_join(res_mr_all_odds_het, plt, by = 'method')
  command3 <- paste0("write.table('mr_results_", wbc, "_", protein, "_all.txt', res_mr_all_odds_het_plt, sep = \"\t\", quote = FALSE, col.names = TRUE, row.names = FALSE)", sep = "")
  eval(parse(text = command3))
  
  ## MR using all pQTLs after Steiger filtering
  if (nrow(filter(dat_all, steiger_dir == "FALSE")) > 0) {
    dat_all_steiger <- filter(dat_all, steiger_dir == "TRUE")
    res_mr_all_steiger <- mr(dat_all_steiger, method_list=c("mr_ivw", "mr_wald_ratio", "mr_egger_regression", "mr_weighted_median"))
    res_mr_all_steiger_odds <- generate_odds_ratios(res_mr_all_steiger)
    if (nrow(filter(dat_all_steiger, remove == "FALSE"))  > 2) {
      het <- mr_heterogeneity(dat_all_steiger, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
    }
    if (nrow(filter(dat_all_steiger, remove == "FALSE"))  == 2) {
      het <- mr_heterogeneity(dat_all_steiger, method_list = c("mr_ivw", "mr_egger_regression"))
    }
    res_mr_all_steiger_odds_het <- left_join(res_mr_all_steiger_odds, het, by = 'method')
    command4 <- paste0("write.table('mr_results_", wbc, "_", protein, "_all_steiger.txt', res_mr_all_steiger_odds_het_plt, sep = \"\t\", quote = FALSE, col.names = TRUE, row.names = FALSE)", sep = "")
    eval(parse(text = command4))
  }
  
  
  ## MR using non_receptor QTLs
  expo_dat_non_receptor <- format_data(cis, type = "exposure", snp_col = "chr_genpos", beta_col = "BETA", se_col = "SE",
                              eaf_col = "EAF", effect_allele_col = "effect_allele", other_allele_col = "other_allele", 
                              pval_col = "log10P", chr_col = "chr", pos_col = "pos", log_pval = TRUE)
  dat_non_receptor <- harmonise_data(exposure_dat = expo_dat_non_receptor, outcome_dat = outco, action = 1)
  res_mr_non_receptor <- mr(dat_non_receptor, method_list=c("mr_ivw", "mr_wald_ratio", "mr_egger_regression", "mr_weighted_median"))
  res_mr_non_receptor_odds <- generate_odds_ratios(res_mr_non_receptor)
  command5 <- paste0("write.table('mr_results_", wbc, "_", protein, "_non_receptor.txt', res_mr_non_receptor_odds, sep = \"\t\", quote = FALSE, col.names = TRUE, row.names = FALSE)", sep = "")
  eval(parse(text = command5))
  
  
  ## MR using receptor QTLs
  expo_dat_recpetor <- format_data(recpetor, type = "exposure", snp_col = "chr_genpos", beta_col = "BETA", se_col = "SE",
                                   eaf_col = "EAF", effect_allele_col = "effect_allele", other_allele_col = "other_allele", 
                                   pval_col = "log10P", chr_col = "chr", pos_col = "pos", log_pval = TRUE)
  dat_recpetor <- harmonise_data(exposure_dat = expo_dat_recpetor, outcome_dat = outco, action = 1)
  res_mr_receptor_each_sum <- NULL
  for (l in 1:nrow(dat_receptor)) {
    command4 <- paste0("dat_receptor_each <- dat_receptor[", l, ",]", sep = "")
    eval(parse(text = command4))
    res_mr_receptor_each <- mr(dat_receptor_each, method_list=c("mr_wald_ratio"))
    res_mr_receptor_each_sum <- rbind(res_mr_receptor_each_sum, res_mr_receptor_each)
  }
  res_mr_receptor_each_odds <- generate_odds_ratios(res_mr_receptor_each_sum)
  command6 <- paste0("write.table('mr_results_", wbc, "_", protein, "_each.txt', res_mr_receptor_each_odds, sep = \"\t\", quote = FALSE, col.names = TRUE, row.names = FALSE)", sep = "")
  eval(parse(text = command6))
}


## Function for small p values
fix_p_function <- function(BETA, SE){
-pchisq((BETA / SE) ^ 2, df = 1, lower.tail = FALSE, log.p = TRUE)/log(10)
}
