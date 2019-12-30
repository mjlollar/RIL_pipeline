#!/bin/bash

##### Align each read

bwa mem dmel_ref_r5.9.fasta $1 $2 | samtools sort -o bam_align.bam

##### Remove Duplicates

samtools rmdup bam_align.bam bam_rmdup.bam

##### Extract first reads in a pair

samtools view -b -f 0X40 bam_rmdup.bam > bam_firsts.bam

###### Combine all processed bams, where *BAM_FIRSTS all bams above

samtools index bam_align.bam
samtools index bam_rmdup.bam
samtools index bam_firsts.bam

#!/bin/bash

samtools index bam_align.bam
samtools index bam_rmdup.bam
samtools index bam_firsts.bam

names="1 2 3"

for name in $names
do
	samtools mpileup -q 20 -r $name bam_align.bam bam_rmdup.bam bam_firsts.bam | gzip - > $name.mpileup.gz
done


##### unzip pileups

gunzip 1.mpileup.gz
gunzip 2.mpileup.gz
gunzip 3.mpileup.gz

##### Extract snps and populate the matrix

##### X Chromosome
perl populate_snp_matrix.pl X.panel < 1.mpileup > X.ahmm_in.panel

##### 2nd Chromosome
chrom2="2L 2R"
for name in $chrom2
do
	perl populate_snp_matrix.pl $name.panel < ${names:1:2}.mpileup > ${names:1:2}.ahmm_in.panel
done

chrom3="3L 3R"
for name in $chrom3
do
	perl populate_snp_matrix.pl $name.panel < ${names:3}.mpileup > ${names:3}.ahmm_in.panel
done

##### Create the ahmm sample file
ls bam_align.bam bam_rmdup.bam bam_firsts.bam | perl -pi -e 's/\n/\t2\n/' > ahmm_in.samples

##### Run AHMM
for name in $names
do
	ancestry_hmm -i $name.ahmm_in.panel -s ahmm_in.samples -a 2 0.5 0.5 -p 0 10000 0.5 -p 1 1 0.5
done