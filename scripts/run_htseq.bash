#!/bin/bash
for sample in DNA
do
     for bam in `ls -1 ~/For_zucker/Blazewicz/Bamfiles/*.bam`
     do
	 htseq-count --format=bam --stranded=no --type=CDS --idattr=ID --mode=union $bam Annotation/AC_18_H184-370_SS-ext-BigMem_MegaHIT_contigs_1k.generated.mapped.gff 1> ${sample}_Counts/`basename ${bam}`.counts &
     done
done

