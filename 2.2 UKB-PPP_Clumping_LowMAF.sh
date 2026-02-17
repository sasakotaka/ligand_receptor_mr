
for i in 19
do

protein=$(cat Protein_Olink_Directory_Size28.txt | awk '{print $1}' | head -"$i" | tail -1)                                        ## protein = CXCL6 (j = 19)

  for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
  do

  n_sig_all=$(cat $protein/"${protein}"_ppp_signif_rsid.txt | awk -v chr=$chr '$1 == chr' | wc -l)

  if [ "$n_sig_all" -gt 0 ]
    then
      plink \
        --bfile $path_reference/$chr \                              ## UKB **   files
        --maf 0.001 \
        --keep ukb.eurFIDIIDPCA.txt \                              ## list of individuals of European ancestry
        --remove w27449_remove.txt \                              ## list of individuals whose informed consent was withdrawn
        --clump $protein/"${protein}"_ppp_signif_rsid.txt \
        --clump-p1 5e-8 \
        --clump-p2 1 \
        --clump-r2 0.001 \
        --clump-kb 1000 \
        --clump-snp-field rsid \
        --clump-field P \
        --memory 30000 \
        --threads 15 \
        --out $protein/"${protein}"_clump_v2_chr"${chr}"

    awk 'NR!=1{print $3}' $protein/"${protein}"_clump_v2_chr"${chr}".clumped >> $protein/"${protein}"_chrall_clump_v2.txt

    fi
  done

  sed -i '/^$/d' $protein/"${protein}"_chrall_clump_v2.txt

done
