import pandas as pd

def output_blocks(filename):
    blocks_list = []
    with open(filename, 'r') as handle:
        for line in handle:
            blocks = line.strip().split('\t')
            if all(block != '.' for block in blocks):
                block_dict = {'id': len(blocks_list), 'blocks': blocks}
                blocks_list.append(block_dict)
                print(block_dict['id'], *blocks)

def output_longest_one_to_two_microsynteny(refbed, qrybed, blocks_filename, output_filename):
    gene_dict = {}
    for bed in [refbed, qrybed]:
        with open(bed, 'r') as handle:
            for line in handle: 
                li = line.strip().split('\t')
                gene = li[3]
                seqid = li[0]
                gene_dict[gene] = seqid
    pd.read_csv(blocks_filename, sep='\t', header=None, names=['ref', 'qry1', 'qry2']) \
        .pipe(lambda df: df[(df['ref'] != '.') & (df['qry2'] != '.')]) \
        .pipe(lambda df: df.assign(
            ref=df['ref'].map(gene_dict),
            qry1=df['qry1'].map(gene_dict),
            qry2=df['qry2'].map(gene_dict)
        )) \
        .value_counts() \
        .to_csv(output_filename, header=['count'])
    # df_blocks = pd.read_csv(blocks_filename, sep='\t', header=None, names=['ref', 'qry1', 'qry2'])
    # df_blocks_filtered = df_blocks[(df_blocks['ref'] != '.') & (df_blocks['qry2'] != '.')]
    # df_seqids = df_blocks_filtered.copy()
    # df_seqids['ref'] = df_blocks_filtered['ref'].map(gene_dict)
    # df_seqids['qry1'] = df_blocks_filtered['qry1'].map(gene_dict)
    # df_seqids['qry2'] = df_blocks_filtered['qry2'].map(gene_dict)
    # df_seqids.value_counts().to_csv(output_filename, header=['count'])

def output_besthit_one_to_two_microsynteny(refbed, qrybed, blocks_filename, anchors_filename, output_filename):
    gene_dict = {}
    for bed in [refbed, qrybed]:
        with open(bed, 'r') as handle:
            for line in handle: 
                li = line.strip().split('\t')
                gene = li[3]
                seqid = li[0]
                gene_dict[gene] = seqid

    ref_blocks = []
    with open(blocks_filename, 'r') as handle:
        for line in handle:
            blocks = line.strip().split('\t')
            if all(block != '.' for block in blocks):
                ref_block = blocks[0]
                ref_blocks.append(ref_block)
    pd.read_csv(anchors_filename, sep='\t', header=None, names=['ref', 'qry', 'bits']) \
        .pipe(lambda df: df[df['ref'].isin(ref_blocks)]) \
        .pipe(lambda df: df.loc[df.groupby('ref')['bits'].idxmax()]) \
        .pipe(lambda df: df.assign(besthit=df['qry'].map(gene_dict))) \
        .to_csv(output_filename, header=['ref', 'qry', 'bits', 'besthit'])

    # df_anchors = pd.read_csv(anchors_filename, sep='\t', header=None, names=['ref', 'qry', 'bits'])
    # df_anchors_filtered = df_anchors[df_anchors['ref'].isin(ref_blocks)]
    # df_anchors_filtered_besthit = df_anchors_filtered.loc[df_anchors_filtered.groupby('ref')['bits'].idxmax()]
    # df_seqids = df_anchors_filtered_besthit.copy()
    # df_seqids['besthit'] = df_anchors_filtered_besthit['qry'].map(gene_dict)
    # df_seqids.to_csv(output_filename, header=['ref', 'qry', 'bits', 'besthit'])



    
