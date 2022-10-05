INPUT=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
REF=/cluster/projects/pughlab/references/TGL/hg38/hg38_random.fa
output=/cluster/projects/pughlab/projects/CHARM/LFS/gatk_coverage/output/swgs
shdir=/cluster/projects/pughlab/projects/CHARM/LFS/gatk_coverage/sh_scripts/swgs
gatk_dir=/cluster/tools/software/gatk/3.8

mkdir -p $shdir
mkdir -p $output

cd $INPUT
ls *.bam > $shdir/bams

cd $shdir
for bam in $(cat bams); do

echo -e "#!/bin/bash\n module load gatk/3.8\n" > $shdir/${bam}.sh
echo -e "java -jar $gatk_dir/GenomeAnalysisTK.jar -T DepthOfCoverage\
 -R $REF\
 -o $output/$bam\
 -I $INPUT/$bam\
 -omitBaseOutput\
 -omitIntervals\
 -omitLocusTable\
 -nt 2\
 -ct 15\
 -pt sample\
 -pt readgroup" >> $shdir/${bam}.sh 
done

cd $shdir

ls *.sh > files
for file in $(cat files);do

sbatch -c 1 -t 24:00:00 --mem 16G $file

done
