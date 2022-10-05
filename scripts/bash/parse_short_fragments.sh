INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG
shdir=/cluster/projects/pughlab/projects/CHARM/LFS/parse_bams/sh_scripts

mkdir -p $INPUTDIR/short_bams
mkdir -p $shdir

cd $INPUTDIR/bams
ls *.bam > bams
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
 cat << EOF > $shdir/${bam}_short.sh
#!/bin/bash
#
module load samtools/1.10

samtools view -h $INPUTDIR/bams/${bam}.bam | \
awk 'substr(\$0,1,1)=="@" || (\$9>= 90 && \$9<=150) || (\$9<=-90 && \$9>=-150)' | \
samtools view -b > $INPUTDIR/short_bams/${bam}_short.bam

cd $INPUTDIR/short_bams

samtools index ${bam}_short.bam

EOF

done 

cd $shdir

ls *_short.sh > files
for file in $(cat files);do
sbatch -c 1 -t 48:00:00 -p all --mem 8G $file

done
