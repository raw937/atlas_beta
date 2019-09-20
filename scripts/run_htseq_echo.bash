echo "#!/bin/bash" > do_htseq
echo "module load python/2.7.8" >> do_htseq
for sample in DNA
do
     for bam in `ls -1 ~/For_zucker/Blazewicz/Bamfiles/*.bam`
     do
	 echo "htseq-count --format=bam --stranded=no --type=CDS --idattr=ID --mode=union $bam Annotation/AC_18_H184-370_SS-ext-BigMem_MegaHIT_contigs_1k.generated.mapped.gff 2> ${sample}_Counts/`basename ${bam}`.err 1> ${sample}_Counts/`basename ${bam}`.counts " >> do_htseq
     done
done
