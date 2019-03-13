import pandas as pd
import numpy as np
import argparse, sys, re
import shlex

cmd = """
for bam in `ls -1 BAM/*.bam`
do
htseq-count --format=bam \
            --stranded=no \
            --type=CDS \
            --idattr=ID \
            --mode=union \
            $bam \
            Annotation/AC_18_H184-370_SS-ext-BigMem_MegaHIT_contigs_1k.mapped.gff \
            > ${bam}.counts &
done 
"""

def map_seqid(bam2gff_file, gff_in, gff_out):
    gff2bam = pd.read_table(bam2gff_file, 
                             names=['gff','bam','length'],index_col='gff')
    with gff_out  as out:
        for line in gff_in:
            m = line.split('\t')
            newline = line
            if len(m) > 8 and m[0] in gff2bam.index:
                newline = re.sub(m[0],  gff2bam.loc[m[0],'bam'], line )
            out.write(newline)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert gff to new mapping')
    parser.add_argument( 'mapfile', nargs='?', 
                         type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument( 'gff_in', 
                         type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('gff_out', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    map_seqid( args.mapfile, args.gff_in, args.gff_out )
