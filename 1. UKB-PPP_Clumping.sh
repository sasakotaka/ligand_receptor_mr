
for i in 2:29
do

protein=$(cat Protein_Olink_Directory_Size28.txt | awk '{print $1}' | head -"$i" | tail -1)
directory=$(cat Protein_Olink_Directory_Size28.txt | awk '{print $3}' | head -"$i" | tail -1)

zcat $directory/*.gz | awk '$13 > 7.3' | awk '$15 = 10^(-$13)' | awk '$16 = $6 * $8' > $protein/"${protein}"_ppp_signif_tentative.txt
## extracting SNPs with genome-wide significance from UKB-PPP summary statistics

cat $protein/"${protein}"_ppp_signif_tentative.txt | awk '{print "chr",$1,":",$2,":",$4,":",$5}' | sed 's/ //g' > $protein/formatID.txt
cat $protein/"${protein}"_ppp_signif_tentative.txt | awk '{print $3}' | sed 's/:/ /g' | awk '{print $2}' > $protein/formatID.txt/position_h37.txt
paste -d " " $protein/"${protein}"_ppp_signif_tentative.txt $protein/formatID.txt $protein/formatID.txt/position_h37.txt > $protein/"${protein}"_ppp_signif.txt

  for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
  do

  n_sig_all=$(cat $protein/"${protein}"_ppp_signif.txt | awk -v chr=$chr '$1 == chr' | wc -l)

  if [ "$n_sig_all" -gt 0 ]
    then
      plink \
        --bfile topmed_chr$chr \                              ## UKB Topmed bgen files
        --maf 0.001 \
        --keep ukb.eurFIDIIDPCA.txt \                              ## list of individuals of European ancestry
        --remove w27449_remove.txt \                              ## list of individuals whose informed consent was withdrawn
        --clump $protein/"${protein}"_ppp_signif.txt \
        --clump-p1 5e-8 \
        --clump-p2 1 \
        --clump-r2 0.001 \
        --clump-kb 1000 \
        --clump-snp-field formatID \
        --clump-field P \
        --memory 30000 \
        --threads 15 \
        --out $protein/"${protein}"_clump_chr"${chr}"

    awk 'NR!=1{print $3}' $protein/"${protein}"_clump_chr"${chr}".clumped >> $protein/"${protein}"_chrall_clump.txt

    fi
  done

  sed -i '/^$/d' $protein/"${protein}"_chrall_clump.txt

done
