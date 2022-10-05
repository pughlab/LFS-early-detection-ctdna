#!/bin/bash

## Extract properly-paired reads for wild-type and mutant seperately
##=== Part I ===
# Link files

projectdir=/cluster/projects/pughlab/projects/CHARM/LFS
dir=/cluster/projects/pughlab/projects/CHARM/LFS/mutation_fragment
INPUT=/cluster/projects/pughlab/external_data/TGL49_CHARM/LFS/LFS_TS/all_unique
mutfile=$dir/CHARM_LFS_mutations_list_formatted.txt
shdir=$dir/sh_scripts
outdir=$dir/output
picard_dir=/cluster/tools/software/picard/2.10.9

mkdir -p $shdir
mkdir -p $outdir

cd $INPUT
ls *bam > $shdir/bam
cd $shdir
sed 's/....$//' bam > bams
rm bam

##=== Part II ===
# Generate bam file including properly-paired reads

for bam in $(cat bams);do
echo $bam
mkdir -p $outdir/$bam/
# Generate .sh
echo -e '#!/bin/bash\n' > $shdir/${bam}.sh
echo -e "module load picard \nmodule load samtools\n" >> $shdir/${bam}.sh
echo -e "cd $outdir/$bam" >> $shdir/${bam}.sh

awk '{print $1}' $mutfile | while read line;do
 gene=$(echo $line | awk -F '.' '{print $1}')
 chr=$(echo $line | awk -F'.' '{print $2}')
 pos=$(echo $line | awk -F'.' '{print $3}')
 ref=$(echo $line | awk -F'.' '{print $4}')
 alt=$(echo $line | awk -F'.' '{print $5}')
 mut=$(echo $line | awk -F'.' '{print $6}')
        
echo -e "\n##=== Make position specific bams ===

##=== ${chr}_${pos}_${ref}_${alt}_${gene}_${mut} ===
        
samtools view -h -f 0x2 $INPUT/${bam}.bam ${chr}:${pos}-${pos} | samtools view -bhSo ${bam}_${gene}_${mut}.bam -
samtools index ${bam}_${gene}_${mut}.bam
samtools view ${bam}_${gene}_${mut}.bam > ${bam}_${gene}_${mut}_reads" >> $shdir/${bam}.sh
        
echo -e "\n##=== Write WT and Mutant bams ===

awk -v pos=$pos -v ref=$ref -v alt=$alt '{loc=pos-\$4+1; split(\$10, arr, \"\"); if(arr[loc]==ref) print \$0}' ${bam}_${gene}_${mut}_reads > ${bam}_${gene}_${mut}_wildtype
awk -v pos=$pos -v ref=$ref -v alt=$alt '{loc=pos-\$4+1; split(\$10, arr, \"\"); if(arr[loc]==alt) print \$0}' ${bam}_${gene}_${mut}_reads > ${bam}_${gene}_${mut}_mutant
uniq ${bam}_${gene}_${mut}_wildtype ${bam}_${gene}_${mut}_wildtypes
uniq ${bam}_${gene}_${mut}_mutant ${bam}_${gene}_${mut}_mutants
mv ${bam}_${gene}_${mut}_wildtypes ${bam}_${gene}_${mut}_wildtype
mv ${bam}_${gene}_${mut}_mutants ${bam}_${gene}_${mut}_mutant" >> $shdir/${bam}.sh
        
echo -e "\n##=== Remove extra files ===

rm ${bam}_${gene}_${mut}_reads
rm ${bam}_${gene}_${mut}.bam*" >> $shdir/${bam}.sh
        
echo -e "\n###=== Next mutation ===" >> $shdir/${bam}.sh
done
        
echo -e "find -type f -empty  -delete" >> $shdir/${bam}.sh
done

# Remove ^M from .sh
cd $shdir
for i in *.sh;do sed -i "s/\r//g" $i;done;

# Check
##=== Part III ===
# Submit scripts

cd $shdir
ls *sh > files
for file in $(cat files); do
sbatch -c 1 -t 2:00:00 -p all --mem 4G $file
done
