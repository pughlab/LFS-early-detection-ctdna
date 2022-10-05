INPUT=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_TS/bams
REF=/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa
output=/cluster/projects/pughlab/projects/CHARM/LFS/gatk_coverage/output/ts_raw
shdir=/cluster/projects/pughlab/projects/CHARM/LFS/gatk_coverage/sh_scripts/ts_raw
interval=/cluster/projects/pughlab/projects/CHARM/intervals/CHARM-hg38_arna.bed

mkdir -p $output
mkdir -p $shdir

cd $INPUT
ls *.bam > $shdir/bams

cd $shdir
for bam in $(cat bams); do

echo -e "#!/bin/bash\n module load gatk/4.1.8.1\n" > $shdir/${bam}.sh
echo -e "gatk DepthOfCoverage\
 -R $REF\
 -O $output/$bam\
 -I $INPUT/$bam\
 -L $interval\
 --omit-interval-statistics false\
 --omit-locus-table true" >> $shdir/${bam}.sh 
done

cd $shdir

ls *.sh > files
for file in $(cat files);do

sbatch -c 1 -t 24:00:00 --mem 8G $file

done
