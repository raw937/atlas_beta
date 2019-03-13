import pandas as pd


def get_counts(counts_file, annotation_file, keys=[['ko_id', 'ko_gene_symbol', 'ko_product', 'ko_ec'],
                                                   ['ko_level1_name'],
                                                   ['ko_level2_name'],
                                                   ['ko_level3_name'],
                                                   ['ko_ec'],
                                                   ['cog_id', 'cog_product'],
                                                   ['cog_level1_name'],
                                                   ['cog_level2_name'],
                                                   ['taxonomy'],
                                                   ['uniprot_ac', 'cazy_id1', 'cazy_class', 'cazy_product', 'cazy_ec', 'cazy_gene_id', 'cazy_taxa', 'cazy_id1', 'cazy_class']]):
    counts = pd.read_csv(counts_file, sep='\t')
    annotation = pd.read_csv(annotation_file, sep='\t')

    merged = counts.merge(annotation, how='left', left_on='gene', right_on='orf')

    out = counts_file.split('.')[0]
    for k in keys:
        tmp = k[:]
        tmp.append('count')
        subset = merged[tmp].copy()
        j = k[0]
        subset.dropna(inplace=True)
        ss_grouped = subset.groupby(k).sum()
        ss_grouped.to_csv(out + '_%s_counts.tsv' % j, sep='\t')


def taxa_counts(counts_file, annotation_file, keys=['ko_id', 'cog_id', 'cazy_id1'], ec=['ko_ec', None, 'cazy_ec']):
    col_names = ['root', 'organism', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    counts = pd.read_csv(counts_file, sep='\t')
    annotation = pd.read_csv(annotation_file, sep='\t')

    out = counts_file.split('.')[0]
    for k, ec in zip(keys, ec):
        print(k)
        merged = counts.merge(annotation, how='left', left_on='gene', right_on='orf')
        if ec is None:
            subset = merged[['taxonomy', k, 'count']].copy()
        else:
            subset = merged[['taxonomy', k, ec, 'count']].copy()
        subset.dropna(inplace=True)
        tmp = subset['taxonomy'].str.split(';', expand=True).ix[:, 0:8]
        tmp.columns = col_names

        for l in [4, 5, 8]:
            print(col_names[l - 1])
            tmp[k] = subset[k]
            if ec is not None:
                tmp[ec] = subset[ec]
            tmp['count'] = subset['count']

            groups = col_names[:l]

            if ec is None:
                groups.append(k)
            else:
                groups.extend([k, ec])
            grouped = tmp.groupby(groups).sum()

            grouped.to_csv(out + '_%s_%s.tsv' % (k, col_names[l - 1]), sep='\t')


if __name__ == '__main__':
    counts_file = '/home/sean/Desktop/21_Unaligned_Extended_Frags_Trimmed_Aligned.bam.counts.CDS.txt'
    annotation_file = '/home/sean/Desktop/21_refseq_eggnog_merged.tsv'
    taxa_counts(counts_file, annotation_file)
