import pandas as pd
import numpy as np
import argparse, sys, re, os
sys.path.append(os.getcwd())
from funtaxacount import get_read_counts

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect read counts into one file')
    parser.add_argument( '--read_count_files', nargs='*', 
                         type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument( '--out', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    read_counts = get_read_counts( args.read_count_files )
    read_counts.to_csv(args.out, sep='\t',index_label='gene')
             
