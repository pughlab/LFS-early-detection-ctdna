INPUTDIR=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_WG/bams
shdir=/cluster/projects/pughlab/projects/CHARM/LFS/insert_size/sh_scripts/swgs
outdir=/cluster/projects/pughlab/projects/CHARM/LFS/insert_size/output/swgs
picard_dir=/cluster/tools/software/picard/2.10.9

mkdir -p $shdir
mkdir -p $outdir

cd $INPUTDIR
ls *.bam > $shdir/bams
cd $shdir
sed 's/....$//' bams > bam
rm bams
mv bam bams

for bam in $(cat bams); do
 cat << EOF > $shdir/${bam}_picard.sh
#!/bin/bash
#
module load picard

java -jar $picard_dir/picard.jar CollectInsertSizeMetrics \
I=$INPUTDIR/${bam}.bam \
O=$outdir/${bam}_picard.txt \
H=$outdir/${bam}.pdf \
M=0 \
W=600 

EOF

done

cd $shdir

ls *picard.sh > files
for file in $(cat files); do
sbatch -c 1 -t 2:00:00 -p all --mem 6G $file
done
