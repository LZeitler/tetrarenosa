## setwd() # make sure wd contains shell scripts invoked here!

require(tidyverse)
require(data.table)

source("ld_easy.R")

prepVCF(
    vcfdir="/nfs/nas22/fs2201/biol_impb_group_evo_euler/data/swiss_tetra/gatk3hc_300_tetra_combined.fi.only4.vcf.gz",
    bcftools.loadpath="module load gcc/4.8.2 gdc perl/5.18.4 samtools/1.9",
    chromosomes=c("scaffold_1","scaffold_2","scaffold_3","scaffold_4","scaffold_5","scaffold_6","scaffold_7","scaffold_8"),
    temp.dir="/cluster/scratch/zeitlerl/pruning"
)

pruneChromosomes(
    temp.dir="/cluster/scratch/zeitlerl/pruning",
    windowsize=10000,
    rlimit=0.1
)

postVCF(
    vcfdir="/nfs/nas22/fs2201/biol_impb_group_evo_euler/data/swiss_tetra/gatk3hc_300_tetra_combined.fi.only4.vcf.gz",
    bcftools.loadpath="module load gcc/4.8.2 gdc perl/5.18.4 samtools/1.9",
    temp.dir="/cluster/scratch/zeitlerl/pruning"
)
