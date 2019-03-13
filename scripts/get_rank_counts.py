import argparse
import pandas as pd
import numpy as np
import os,sys, re
sys.path.append(os.getcwd())
from funtaxacount import merge_orf_and_funtax, get_rank_counts

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate rank counts')
    parser.add_argument( 'funtax', nargs='?',type=argparse.FileType('r'))
    parser.add_argument( 'orf', type=argparse.FileType('r'))
    parser.add_argument( 'read_counts', type=argparse.FileType('r'))
    parser.add_argument( 'ncbi_megan_map', type=argparse.FileType('r') )
    parser.add_argument( 'ncbi_tree', type=argparse.FileType('r') )
    parser.add_argument( 'merge_map', type=argparse.FileType('r') )
    parser.add_argument( '--rank', choices=['phylum','class'], default='phylum') 
    parser.add_argument( '--out', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    args = parser.parse_args()
    rank_counts = get_rank_counts( args.funtax.name, args.orf.name,  args.read_counts, args.ncbi_tree, args.ncbi_megan_map,  args.merge_map, args.rank )
    rank_counts.to_csv(args.out, sep='\t')
        
    
