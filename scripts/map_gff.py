import pandas as pd
import numpy as np
import argparse, sys, re
from funtaxacount import map_seqid

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert gff to new mapping')
    parser.add_argument( 'mapfile', nargs='?', 
                         type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument( 'gff_in', 
                         type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument( 'gff_out', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    map_seqid( args.mapfile, args.gff_in, args.gff_out )
