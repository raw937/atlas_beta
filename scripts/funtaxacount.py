import pandas as pd
import numpy as np
import argparse, sys, re

orf_names = ['ORF_ID', 'Contig', 'COG', 'KO'] #, 'product']

def merge_orf_and_funtax( orf_file, funtax_file ):
    """
    Takes an orf file and a funtaxa file and returns the merge
    """
    orf_df = pd.read_table(orf_file, header=None, names=orf_names, index_col='ORF_ID', usecols=orf_names, engine='python', encoding="ISO-8859-1", quoting=3)
    funtax_df = pd.read_table(funtax_file, index_col='ORF_ID', engine='python', encoding="ISO-8859-1", quoting=3)
    funtax_df[['COG','KO']] = orf_df[['COG','KO']]
    funtax_df['taxonId'] = funtax_df['taxonomy'].replace(r'.+\(([0-9]+)\)', value=r'\1', regex=True)
    genes = funtax_df.reset_index()
    genes['gene'] = genes['ORF_ID']
    return genes.set_index('gene')

def get_function_counts( function, orf_file, funtax_file, read_count_file ):
    genes = merge_orf_and_funtax( orf_file, funtax_file )
    read_counts = pd.read_table( read_count_file, index_col='gene', engine='python')
    read_counts[function] = genes[function]
    function_counts = read_counts.groupby(function).sum()
    #if function = 'ec':
    #    ecs_df = pd.read_table(ec_file, index_col='EC')
    #    return pd.concat([ec_df, function_counts], axis=1, join_axes=[function_counts.index])
    return function_counts


def reindex_contig_gene_id( contig_id ):
    return lambda match_obj: 'ID={}_{};'.format( contig_id, int(match_obj.group(2)) - 1 )

def map_seqid(annotation2assembly_file, gff_in, gff_out, map_id_p = False):
    cgRE = re.compile(r'ID=([0-9]+)_([0-9]+);')    
    annotation2assembly_map = pd.read_table(annotation2assembly_file,
                            names=['annotation','assembly','length'],index_col='annotation')
    with gff_out  as out:
        for line in open(gff_in.name, 'r', encoding='ISO-8859-1'):
            m = line.split('\t')
            newline = line
            if len(m) > 8 and m[0] in annotation2assembly_map.index:
                newline = annotation2assembly_map.loc[m[0],'assembly'] + '\t' +  '\t'.join(m[1:])
                if map_id_p:
                    newline = cgRE.sub( reindex_contig_gene_id( m[0] ), newline )
            out.write(newline)

            
def generate_gff( mapfile, funtax_orf_file ):
    """
    Takes the mapfile and annotation file and generates a gff file that maps reads in the bamfile to genes
    """
    annotation2assembly_map = pd.read_table(mapfile,
                                            names=['annotation','assembly','length'],
                                            index_col='annotation')
    funtax_gff = pd.read_table( funtax_orf_file.name, engine='python', encoding='ISO-8859-1', quoting=3)
    funtax_gff['seqid'] = funtax_gff.join(annotation2assembly_map, on='Contig_Name')['assembly']
    funtax_gff['source'] = 'Prodigal_v2.00'
    funtax_gff['type'] =  'CDS'
    funtax_gff['score'] = 100.0
    funtax_gff['phase'] =  0
    funtax_gff['attributes'] = funtax_gff['ORF_ID'].str.replace(r'(.*)', r'ID=\1;')
    return funtax_gff[['seqid','source', 'type','start', 'end', 'score', 'strand','phase','attributes']]


def get_funtaxa_counts( funtax_file, orf_file, read_count_file, ncbi_tree_file, ncbi_megan_map_file, merged_file, rank, function ):
    ncbi_tree = get_ncbi_tree( ncbi_tree_file )
    ncbi_megan_map = get_ncbi_megan_map( ncbi_megan_map_file )
    genes = merge_orf_and_funtax( orf_file, funtax_file )
    read_counts = pd.read_table( read_count_file, index_col='gene', engine='python')
    read_counts['taxonId'] = check_merged( genes['taxonId'], merged_file )
    read_counts[function] = genes[function]
    fun_taxon_counts = read_counts.groupby(['taxonId', function]).sum().reset_index()
    fun_taxon_counts[rank] = fun_taxon_counts['taxonId'].apply( lambda x: get_rank( ncbi_tree, x, rank ) )
    fun_rank_counts = fun_taxon_counts.drop('taxonId', axis=1).groupby([rank, function]).sum()
    fun_rank_counts.insert(0, 'ancestry', pd.Series([get_ancestry(ncbi_tree, r, ncbi_megan_map )
                                                     for (r, f) in fun_rank_counts.index],
                                                    index=fun_rank_counts.index))
    return fun_rank_counts.dropna()

def get_read_counts(read_count_files ):
    read_count_array = []
    header = []
    for read_count_file in read_count_files:
        header.append(read_count_file.name)
        read_count_array.append(pd.read_table( read_count_file, index_col=0,header=None, names=['gene',header[-1]]))
    return pd.concat(read_count_array, axis=1)

def check_merged( taxonId, merged_file ):
    merged_dict = get_merged_map( merged_file )
    taxonId = pd.to_numeric(taxonId, errors= 'coerce')
    return taxonId.apply( lambda x: merged_dict.get(x,x))

def get_rank_counts( funtax_file, orf_file, read_count_file, ncbi_tree_file, ncbi_megan_map_file, merged_file, rank ):
    ncbi_tree = get_ncbi_tree( ncbi_tree_file )
    ncbi_megan_map = get_ncbi_megan_map( ncbi_megan_map_file )
    genes = merge_orf_and_funtax( orf_file, funtax_file )
    read_counts = pd.read_table( read_count_file, index_col='gene', engine='python')
    read_counts['taxonId'] = check_merged( genes['taxonId'], merged_file )
    taxon_counts = read_counts.groupby('taxonId').sum().reset_index()
    taxon_counts[rank] = taxon_counts['taxonId'].apply( lambda x: get_rank( ncbi_tree, x, rank ) )
    rank_counts = taxon_counts.drop('taxonId', axis=1).groupby( rank ).sum().reset_index()
    #Insert ancestry into this column
    rank_counts.insert(0, 'ancestry', rank_counts[rank].apply( lambda x: get_ancestry(ncbi_tree, x, ncbi_megan_map ) ) )
    return rank_counts.set_index(rank).dropna()

def get_merged_map( merged_file ):
    return  pd.read_table( merged_file ,
                     header=None, sep=r'\t\|',
                     names=['tax_id','merged_tax_id'], 
                     engine='python',
                     index_col=0,usecols=[0,1]).to_dict()['merged_tax_id']

def get_ncbi_megan_map( meganfile ):
    ncbi_megan_map = {}
    for line in meganfile:
        fields = line.split("\t")
        taxonId, taxonName = int(fields[0].strip()), fields[1].strip()
        ncbi_megan_map[taxonId] = taxonName
    return ncbi_megan_map

def get_ncbi_tree( ncbi_tree_file ):
    return  pd.read_table( ncbi_tree_file,
                     header=None, sep=r'\t\|\t',
                     names=['tax_id','parent tax_id','rank'], 
                     engine='python',
                     index_col=0,usecols=[0,1,2])
def get_rank( ncbi_tree, tax_id, rank):
    old_id = -1
    try:
        current_id = int(tax_id)
    except ValueError:
        return tax_id
    while current_id in ncbi_tree.index and current_id != old_id:

        if ncbi_tree.loc[current_id, 'rank'] == rank:
            return current_id
        else:
            old_id = current_id
            current_id = ncbi_tree.loc[current_id, 'parent tax_id']
    return current_id

def get_ancestry( ncbi_tree, taxon_id, ncbi_megan_map ):
    ancestry = []
    old_id = -1
    try:
        current_id = int(taxon_id)
    except ValueError:
        return np.nan
    while current_id in ncbi_tree.index and current_id != old_id:
        old_id = current_id
        current_id = ncbi_tree.loc[current_id, 'parent tax_id']
        ancestry.insert(0, old_id )
    return ';'.join([ ncbi_megan_map[i] for i in ancestry if i in ncbi_megan_map])
        
