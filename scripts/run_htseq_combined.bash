#!/bin/bash
for sample in Combined_RNA_BAM Combined_DNA_BAM
do
     for bam in `ls -1 ${sample}_BAM/*no-glycine*.bam`
     do
	 htseq-count --format=bam --stranded=no --type=CDS --idattr=ID --mode=union $bam Assembly/White-moleculo_LongRead.generated.mapped.gff 1> ${sample}_Counts/`basename ${bam}`.counts &
     done
done

