PLINK=/media/NeoScreen/Software/stow/plink_1.07/bin/plink
GTOOL=/media/NeoScreen/Software/gtool

#INPUT is the full path + file root to the .ped and .map files split by chr
INPUT=/media/NeoScreen/NeSc_home/CHSH/Trios/POE/miRNA_PoE/INPUT/MERGED_chr22

ANNO=/media/NeoScreen/NeSc_home/CHSH/Trios/POE/miRNA_PoE/annotation.RData

# Expression data set
EXP=/media/NeoScreen/NeSc_home/CHSH/Trios/POE/miRNA_PoE/expression.RData

# OUTPUT filepath (without / at the end)
OUTPUT=/media/NeoScreen/NeSc_home/CHSH/Trios/POE/miRNA_PoE/OUTPUT


SHAPEIT=/media/NeoScreen/Software/bin/shapeit

-M /media/NeoScreen/NeSc_home/CHSH/Trios/POE/miRNA_PoE/genetic_map_b37/genetic_map_chr${CHR}_combined_b37.txt