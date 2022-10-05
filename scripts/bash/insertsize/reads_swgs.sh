INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/sWGS_bams_2021-05-12
shdir=/cluster/projects/pughlab/projects/CHARM/LFS/insert_size/sh_scripts/swgs
outdir=/cluster/projects/pughlab/projects/CHARM/LFS/insert_size/output/swgs

mkdir -p $shdir
mkdir -p $outdir

cd $INPUTDIR
ls *.bam > bams
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
 cat << EOF > $shdir/${bam}_reads.sh
#!/bin/bash
#
module load samtools

samtools view -c -F 260 $INPUTDIR/${bam}.bam > $outdir/${bam}_reads

EOF

done

cd $shdir

ls *reads.sh > files
for file in $(cat files); do
sbatch -c 1 -t 2:00:00 -p all --mem 4G $file
done
