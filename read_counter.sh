#!/bin/bash

#### Ensure you are in the directory that contains your fastq files
#### Outputs space-delinated text file containing read count, fastq filename

for file in *.fastq.gz
do
  ##### For decompressed files use a=$(echo $(cat $file|wc -l)/4|bc)
  a=$(echo $(zcat < $file|wc -l)/4|bc)
	echo $a $file >> readcounts.txt
done 
