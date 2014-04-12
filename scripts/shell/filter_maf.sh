#!/bin/sh
# filter maf file from "cancer gene saturation" paper
cat PanCan.maf | awk -F"\t" '$5!="IGR" && $5!="Unknown" && $5!="Intron" && $5!="RNA" && $5!="Non-coding_Transcript" && $3!="Unknown" && $5!~/^[35]'\''/' > PanCan.filtered.maf
# conver filtered maf to cravat format
python scripts/python/maf2cravat.py PanCan.filtered.maf PanCan.filtered.cravat.txt

# filter maf file from "haploinsufficiency and triplosensitivity" paper
cat Mutation_Dataset.txt | awk -F"\t" '$6!~/^[35]'\''/ && $6!="Intron"' > Mutation_Dataset.filtered.txt 
