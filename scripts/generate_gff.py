import pandas as pd
import numpy as np
import argparse, sys, re
from funtaxacount import generate_gff

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert gff to new mapping')
    parser.add_argument( 'mapfile', nargs='?', 
                         type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument( 'annotation',
                         type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument( 'gff_out', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    generate_gff( args.mapfile, args.annotation ).to_csv(args.gff_out, sep='\t', index=False, header=False )
