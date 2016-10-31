#!/bin/bash -v
#$ -cwd -V
#$ -j y
#$ -P lindgren.prjc -q short.qc
#$ -t 22-1
##(no limit on jobs) -tc 150

echo "************************************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "************************************************************"


SHAPEIT=/apps/well/shapeit/2.r790/shapeit

## for running the script on a head node, provide index as first argument
if [ -z "$SGE_TASK_ID" ]; then SGE_TASK_ID=$1; fi

CHR=$SGE_TASK_ID

$SHAPEIT --input-ped /well/lindgren/alexd/HAPMAP/hapmap3_pop/hapmap3_r2_updateb37_fwd.MERGED.chr${CHR}.ped /well/lindgren/alexd/HAPMAP/hapmap3_pop/hapmap3_r2_updateb37_fwd.MERGED.chr${CHR}.map \
        -M /well/lindgren/alexd/1000GP_Phase3/genetic_map_chr${CHR}_combined_b37.txt \
        --duohmm \
        -W 5 \
        -O /well/lindgren/alexd/HAPMAP/SHAPEIT/MERGED_chr${CHR}
        --exclude-snp /well/lindgren/alexd/HAPMAP/hapmap3_pop/hapmap3_r2_updateb37_fwd.MERGED.chr${CHR}.dupvar
        #--output-max  \
        #--output-graph /panfs/panasas01/sscm/xt610/TRIOS/IMPUTE/SHAPEIT/TRIO-graph


## "This produces duohmm-example-corrected.haps which is a standard SHAPEIT2 format haplotype file that has been corrected for pedigree structure. An additional benefit is that a child's haplotypes are now ordered with the paternal haplotype first and the maternal haplotype second."
