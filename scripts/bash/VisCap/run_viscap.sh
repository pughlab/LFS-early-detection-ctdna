#!/bin/bash

#SBATCH --mem=4G

module load R/3.6.1

Rscript /cluster/projects/pughlab/projects/CHARM/LFS/scripts/VisCap/VisCap-arna.R /cluster/projects/pughlab/projects/CHARM/LFS/gatk_coverage/output/all_unique  /cluster/projects/pughlab/projects/CHARM/LFS/VisCap /cluster/projects/pughlab/projects/CHARM/LFS/scripts/VisCap.cfg
