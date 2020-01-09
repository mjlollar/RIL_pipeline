#!/bin/bash


for RIL in $(cat RIL_list.txt)
do
	##### Unzip fastq file
	gunzip ${RIL}_R1.fastq.gz
	gunzip ${RIL}_R2.fastq.gz

	##### Align each read
	bwa mem dmel_ref_r5.9.fasta ${RIL}_R1.fastq ${RIL}_R2.fastq | /mnt/sas0/AD/mlollar/samtools sort -o bam_align.bam
	
	##### Remove Duplicates
	/mnt/sas0/AD/mlollar/samtools rmdup bam_align.bam bam_rmdup.bam
	
	##### Extract first reads in a pair
	/mnt/sas0/AD/mlollar/samtools view -b -f 0X40 bam_rmdup.bam > bam_firsts.bam
	
	###### Combine all processed bams, where *BAM_FIRSTS all bams above
	/mnt/sas0/AD/mlollar/samtools index bam_align.bam
	/mnt/sas0/AD/mlollar/samtools index bam_rmdup.bam
	/mnt/sas0/AD/mlollar/samtools index bam_firsts.bam
	
	names="1 2 3"
	
	for name in $names
	do
		/mnt/sas0/AD/mlollar/samtools mpileup -q 20 -r $name bam_align.bam bam_rmdup.bam bam_firsts.bam | gzip - > ${name}.mpileup.gz
	done


	##### unzip pileups
	gunzip 1.mpileup.gz
	gunzip 2.mpileup.gz
	gunzip 3.mpileup.gz

	##### Extract snps and populate the matrix
	##### X Chromosome
	perl mnt/sas0/AD/mlollar/bin/populate_snp_matrix.pl X.panel < 1.mpileup > X.ahmm_in.panel
	
	##### 2nd Chromosome
	chrom2="2L 2R"
	for name in $chrom2
	do
		perl mnt/sas0/AD/mlollar/bin/populate_snp_matrix.pl ${name}.panel < 2.mpileup > ${name}.ahmm_in.panel
	done
	
	##### 3rd Chromosome
	chrom3="3L 3R"
	for name in $chrom3
	do
		perl mnt/sas0/AD/mlollar/bin/populate_snp_matrix.pl ${name}.panel < 3.mpileup > ${name}.ahmm_in.panel
	done
	
	##### Create the ahmm sample file
	ls bam_align.bam bam_rmdup.bam bam_firsts.bam | perl -pi -e 's/\n/\t2\n/' > ahmm_in.samples

	chroms="X 2L 2R 3L 3R"

	##### Run AHMM
	for name in $chroms
	do
		ancestry_hmm -i ${name}.ahmm_in.panel -s ahmm_in.samples -a 2 0.5 0.5 -p 0 10000 0.5 -p 1 1 0.5
		mv bam_rmdup.bam.posterior ../ahmm_outputs/${RIL}.${name}.bam_rmdup.bam.posterior
		mv bam_firsts.bam.posterior ../ahmm_outputs/${RIL}.${name}.bam_firsts.bam.posterior
		mv bam_align.bam.posterior ../ahmm_outputs/${RIL}.${name}.bam_align.bam.posterior
	done

	##### Remove unzipped file to conserve disk space
	rm ${RIL}_R1.fastq
	rm ${RIL}_R2.fastq

done

