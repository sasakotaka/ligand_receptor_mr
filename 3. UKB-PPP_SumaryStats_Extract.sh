
for i in 2:29
do

protein=$(cat Protein_Olink_Directory_Size28.txt | awk '{print $1}' | head -"$i" | tail -1)
genpos_list=$(cat $protein/"${protein}"_chrall_clump.txt | awk '$2 = "-e"' | awk '{print $2,$1}' ORS=' ')

echo -e "cat $protein/"$protein"_ppp_signif.txt | grep -w "$genpos_list" > $protein/"$protein"_ppp_signif_extract.txt" > newfile_"$protein"
chmod u+x newfile_"$protein"
./newfile_"$protein"

done


for i in 19
do

protein=$(cat Protein_Olink_Directory_Size28.txt | awk '{print $1}' | head -"$i" | tail -1)
rsid_list=$(cat $protein/"${protein}"_chrall_clump_v2.txt | awk '$2 = "-e"' | awk '{print $2,$1}' ORS=' ')

echo -e "cat $protein/"$protein"_ppp_signif_rsid.txt | grep -w "$rsid_list" > $protein/"$protein"_ppp_signif_extract.txt" > newfile_"$protein"_v2
chmod u+x newfile_"$protein"_v2
./newfile_"$protein"_v2

done
