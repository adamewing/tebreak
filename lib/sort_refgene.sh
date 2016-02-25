zcat refGene.txt.gz | sort -k3,3 -k5,5n > refGene_sorted.txt
bgzip refGene_sorted.txt
tabix -s 3 -b 5 -e 6 refGene_sorted.txt.gz
