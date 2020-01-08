# RIL_pipeline

Shell script to: 
 - Align Recombinant Inbred Line sequences to D. melanogaster reference (v5.9). Illumina Nova-Seq (2x150bp) paired-end reads.
 - Process bams, obtain SNP information, generate ahmm files
 - Perform ancestry inferences using Ancestry_HMM[1].

Script requirements:
 - Genome version 5.9.
 - samtools version 1.6
 - bwa version 0.7.17-r1188
 - Ancestry_HMM version 0.94
 - Perl script populate_SNP_matrix.pl

References:

[1] Corbett-Detig, R. and Nielsen, R., 2017. A hidden Markov model approach for simultaneously estimating
local ancestry and admixture time using next generation sequence data in samples of arbitrary
ploidy. PLoS genetics, 13(1), p.e1006529.
