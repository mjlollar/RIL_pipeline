#!/bin/bash


for RIL in $(cat RIL_list.txt)
do
	##### Unzip fastq file
	zcat ${RIL}_R1.fastq.gz > ${RIL}_R1.fastq
	zcat ${RIL}_R2.fastq.gz > ${RIL}_R2.fastq

	##### Align each read
	bwa mem dmel_ref_r5.9.fasta ${RIL}_R1.fastq ${RIL}_R2.fastq | /mnt/sas0/AD/mlollar/bin/samtools-1.6/samtools sort -o bam_align.bam
	
	##### Remove Duplicates
	/mnt/sas0/AD/mlollar/bin/samtools-1.6/samtools rmdup bam_align.bam bam_rmdup.bam
	
	##### Extract first reads in a pair
	/mnt/sas0/AD/mlollar/bin/samtools-1.6/samtools view -b -f 0X40 bam_rmdup.bam > bam_firsts.bam
	
	###### Combine all processed bams, where *BAM_FIRSTS all bams above
	/mnt/sas0/AD/mlollar/bin/samtools-1.6/samtools index bam_align.bam
	/mnt/sas0/AD/mlollar/bin/samtools-1.6/samtools index bam_rmdup.bam
	/mnt/sas0/AD/mlollar/bin/samtools-1.6/samtools index bam_firsts.bam
	
	# Each number in names corresponds to chromosome region, respectively: X, 2L, 2R, 3L, 3R
	names="4 3 7 5 8"
	
	for name in $names
	do
		/mnt/sas0/AD/mlollar/bin/samtools-1.6/samtools mpileup -q 20 -r $name bam_align.bam bam_rmdup.bam bam_firsts.bam | gzip - > ${name}.mpileup.gz
	done

	##### Extract snps and populate the matrix
	##### X Chromosome
	gunzip 4.mpileup.gz
	perl populate_snp_matrix.pl X.panel < 4.mpileup > X.ahmm_in.panel
	
	##### 2L
	gunzip 3.mpileup.gz
	perl populate_snp_matrix.pl 2L.panel < 3.mpileup > 2L.ahmm_in.panel
	
	##### 2R
	gunzip 7.mpileup.gz
	perl populate_snp_matrix.pl 2R.panel < 7.mpileup > 2R.ahmm_in.panel
	
	##### 3L
	gunzip 5.mpileup.gz
	perl populate_snp_matrix.pl 3L.panel < 5.mpileup > 3L.ahmm_in.panel

	##### 3R
	gunzip 8.mpileup.gz
	perl populate_snp_matrix.pl 3R.panel < 8.mpileup.gz > 3R.ahmm_in.panel
	
	##### Create the ahmm sample file
	ls bam_align.bam bam_rmdup.bam bam_firsts.bam | perl -pi -e 's/\n/\t2\n/' > ahmm_in.samples

	chroms="X 2L 2R 3L 3R"

	##### Run AHMM
	for name in $chroms
	do
		/bin/Ancestry_HMM/src/ancestry_hmm -i ${name}.ahmm_in.panel -s ahmm_in.samples -a 2 0.5 0.5 -p 0 10000 0.5 -p 1 1 0.5
		mv bam_rmdup.bam.posterior ../ahmm_outputs/${RIL}.${name}.bam_rmdup.bam.posterior
		mv bam_firsts.bam.posterior ../ahmm_outputs/${RIL}.${name}.bam_firsts.bam.posterior
		mv bam_align.bam.posterior ../ahmm_outputs/${RIL}.${name}.bam_align.bam.posterior
	done

	##### Remove unzipped file to conserve disk space
	rm ${RIL}_R1.fastq
	rm ${RIL}_R2.fastq
	##### Remove Pileups to avoid overwrite prompt
	rm 3.mpileup
	rm 4.mpileup
	rm 5.mpileup
	rm 7.mpileup
	rm 8.mpileup
done
