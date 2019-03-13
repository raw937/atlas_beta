import pandas as pd
import os
import argparse


def get_counts(counts_file, annotation_file, keys=[['ko_id', 'ko_gene_symbol', 'ko_product', 'ko_ec'],
                                                   ['ko_level1_name'],
                                                   ['ko_level2_name'],
                                                   ['ko_level3_name'],
                                                   ['ko_ec'],
                                                   ['cog_id', 'cog_product'],
                                                   ['cog_level1_name'],
                                                   ['cog_level2_name'],
                                                   ['taxonomy'],
                                                   ['uniprot_ac', 'cazy_id1', 'cazy_product', 'cazy_ec', 'cazy_gene_id', 'cazy_taxa', 'cazy_class'],
                                                   ['taxonomy', 'refseq_product'],
                                                   ['taxonomy', 'ko_id'],
                                                   ['taxonomy', 'ko_ec'],
                                                   ['taxonomy', 'cog_id'],
                                                   ['taxonomy', 'cazy_id1'],
                                                   ['cazy_id1']]):


    print("Reading counts...")
    counts = pd.read_csv(counts_file, sep='\t', skiprows=1, dtype=str)

    print("Reading annotations...")
    annotation = pd.read_csv(annotation_file, sep='\t', dtype=str)

    print("Preprocessing...")
    counts = counts.drop(counts.columns[[1, 2, 3, 4, 5]], axis=1)
    counts.ix[:, 1:] = counts.ix[:, 1:].astype(float)

    merged = counts.merge(annotation, how='left', left_on='Geneid', right_on='orf')

    samples = counts.columns[1:]

    out = os.path.splitext(counts_file)[0]
    c = len(keys)
    for i, k in enumerate(keys):
        print("Generating file %d/%d..." % (i + 1, c))
        tmp = k[:]
        tmp.extend(samples)
        subset = merged[tmp].copy()
        j = '-'.join(k)
        subset.dropna(inplace=True)

        ss_grouped = subset.groupby(k).sum()
        ss_grouped.to_csv(out + '_%s.tsv' % j, sep='\t')
    print("Done.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parses quantitation files.')

    # required
    parser.add_argument('counts', help='Path to counts file.')
    parser.add_argument('annotation', help='Path to annotation file.')

    args = parser.parse_args()

    get_counts(args.counts, args.annotation)
