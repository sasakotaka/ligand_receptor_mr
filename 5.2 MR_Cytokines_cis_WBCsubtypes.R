
## Load libraries

library('vroom')
library('dplyr')
library('tidyverse')
library('glue')
library('TwoSampleMR')


## Load data

file <- 'Protein_Olink_Directory_Size_Chr28.txt'
protein_list28 <- vroom(file, col_name = TRUE)

file <- 'outco_dat_neutrophil_formatted.txt'
outco_neutro <- vroom(file, col_names = TRUE)

file <- 'outco_dat_monocyte_formatted.txt'
outco_mono <- vroom(file, col_names = TRUE)

file <- 'outco_dat_lymphocyte_formatted.txt'
outco_lymph <- vroom(file, col_names = TRUE)

wbc_list <- c("neutro", "mono", "lymph")


## Perform MR

for (i in 1:3){
  wbc = wbc_list[i]
  command1 <- paste0("outco <- outco_", wbc, sep = "")
  eval(parse(text = command1))
  dat_cis_rbind <- NULL
  res_mr_cis_rbind <- NULL

  for (j in 1:27) {
    protein = protein_list28$Protein[j]
    protein_chr = protein_list28$Chr[j]
    protein_tss = protein_list28$TSS[j]
    
    command2 <- paste0("file <- '", protein, "_ppp_signif_extract.txt'", sep = "")
    eval(parse(text = command2))
    
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
    cis <- filter(all, CIS_TRANS == "cis")
    
    expo_dat_cis <- format_data(cis, type = "exposure", snp_col = "chr_genpos", beta_col = "BETA", se_col = "SE",
                                eaf_col = "EAF", effect_allele_col = "effect_allele", other_allele_col = "other_allele", 
                                pval_col = "log10P", chr_col = "chr", pos_col = "pos", log_pval = TRUE)
    dat_cis <- harmonise_data(exposure_dat = expo_dat_cis, outcome_dat = outco, action = 1)
    dat_cis_rbind <- rbind(dat_cis_rbind, dat_cis)
    res_mr_cis <- mr(dat_cis, method_list=c("mr_ivw", "mr_wald_ratio", "mr_egger_regression", "mr_weighted_median"))
    res_mr_cis_odds <- generate_odds_ratios(res_mr_cis)
    res_mr_cis_rbind <- rbind(res_mr_cis_rbind, res_mr_cis_odds)
  }
  res_mr_cis_rbind$bonf_P <- p.adjust(res_mr_cis_rbind$pval, method = "bonferroni")
  command3 <- paste0("write.table('mr_data_", wbc, "_cytokines_cis.txt', dat_cis_rbind, sep = \"\t\", quote = FALSE, col.names = TRUE, row.names = FALSE)
                     write.table('mr_results_", wbc, "_cytokines_cis.txt', res_mr_cis_rbind, sep = \"\t\", quote = FALSE, col.names = TRUE, row.names = FALSE)", sep = "")
  eval(parse(text = command3))
}
