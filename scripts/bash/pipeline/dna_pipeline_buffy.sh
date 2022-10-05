#!/bin/bash

module load perl
module load samtools/1.9

perl /cluster/projects/pughlab/projects/CHARM/pipeline-suite/pughlab_dnaseq_pipeline.pl \
-t /cluster/projects/pughlab/projects/CHARM/LFS/configs/LFS_pipeline_bc.yaml \
-d /cluster/projects/pughlab/projects/CHARM/LFS/configs/LFS_bams_bc.yaml \
--variant_calling \
-c slurm \
--remove
