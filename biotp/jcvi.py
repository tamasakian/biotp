import pandas as pd

def output_blocks(filename):
    blocks_list = []
    with open(filename, 'r') as handle:
        for line in handle:
            blocks = line.strip().split('\t')
            if all(block != '.' for block in blocks):
                block_dict = {'id': len(blocks_list) + 1, 'blocks': blocks}
                blocks_list.append(block_dict)
                print(block_dict['id'], *blocks)

def output_longest_one_to_one_microsynteny(refbed, qrybed, blocks_filename, output_filename):
    gene_dict = {}
    for bed in [refbed, qrybed]:
        with open(bed, 'r') as handle:
            for line in handle: 
                li = line.strip().split('\t')
                gene = li[3]
                seqid = li[0]
                gene_dict[gene] = seqid
    pd.read_csv(blocks_filename, sep='\t', header=None, names=['ref', 'qry1']) \
        .pipe(lambda df: df[(df['ref'] != '.') & (df['qry1'] != '.')]) \
        .pipe(lambda df: df.assign(
            ref=df['ref'].map(gene_dict),
            qry1=df['qry1'].map(gene_dict))) \
        .value_counts() \
        .to_csv(output_filename, header=['count'])

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
        .pipe(lambda df: df[(df['ref'] != '.') & (df['qry1'] != '.') & (df['qry2'] != '.')]) \
        .pipe(lambda df: df.assign(
            ref=df['ref'].map(gene_dict),
            qry1=df['qry1'].map(gene_dict),
            qry2=df['qry2'].map(gene_dict))) \
        .value_counts() \
        .to_csv(output_filename, header=['count'])

def output_longest_one_to_three_microsynteny(refbed, qrybed, blocks_filename, output_filename):
    gene_dict = {}
    for bed in [refbed, qrybed]:
        with open(bed, 'r') as handle:
            for line in handle: 
                li = line.strip().split('\t')
                gene = li[3]
                seqid = li[0]
                gene_dict[gene] = seqid
    pd.read_csv(blocks_filename, sep='\t', header=None, names=['ref', 'qry1', 'qry2', 'qry3']) \
        .pipe(lambda df: df[(df['ref'] != '.') & (df['qry1'] != '.') & (df['qry2'] != '.') & (df['qry3'] != '.')]) \
        .pipe(lambda df: df.assign(
            ref=df['ref'].map(gene_dict),
            qry1=df['qry1'].map(gene_dict),
            qry2=df['qry2'].map(gene_dict),
            qry3=df['qry3'].map(gene_dict))) \
        .value_counts() \
        .to_csv(output_filename, header=['count'])

def output_besthit_one_to_two_microsynteny(refbed, qrybed, blocks_filename, anchors_filename, output_filename):
    gene_dict = {}
    for bed in [refbed, qrybed]:
        with open(bed, 'r') as handle:
            for line in handle: 
                li = line.strip().split('\t')
                if len(li) != 6:
                    continue
                seqid, start, end, name, score, strand = li
                gene_dict[name] = seqid

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

def output_besthit_one_to_two_synteny(refbed, qrybed, blocks_filename, anchors_filename, output_filename):
    gene_dict = {}
    for bed in [refbed, qrybed]:
        with open(bed, 'r') as handle:
            for line in handle: 
                li = line.strip().split('\t')
                if len(li) != 6:
                    continue
                seqid, start, end, name, score, strand = li
                gene_dict[name] = seqid

    refs = []
    with open(blocks_filename, 'r') as handle:
        for line in handle:
            li = line.strip().split('\t')
            if len(li) != 3:
                continue
            ref, qry1, qry2 = li
            if ref == "." or qry1 == "." or qry2 == ".":
                continue
            refs.append(ref)

    hits = []
    with open(anchors_filename, 'r') as handle:
        for line in handle:
            li = line.strip().split('\t')
            if len(li) != 3:
                continue
            ref, qry, bits = li
            bits = float(bits)
            if ref not in refs:
                continue
            hits.append((ref, qry, bits))

    with open(output_filename, 'w') as outfile:
        outfile.write("ref_seq\tref\tqry_seq\tqry\tbits\n") 
        for ref, qry, bits in hits:
            ref_seq = gene_dict[ref]
            qry_seq = gene_dict[qry]
            outfile.write(f"{ref_seq}\t{ref}\t{qry_seq}\t{qry}\t{bits}\n")



    
