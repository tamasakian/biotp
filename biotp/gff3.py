#!/usr/bin/env python3

"""Library for processing GFF3 files.

Functions
---------
output_seqid_strand_locs_by_gene_id: Output seqid, strand, start, and end by partial match of gene id.

"""


def output_seqid_strand_locs_by_gene_id(filename, gene_id):
    """Output seqid, strand, start, and end by partial match of gene id.

    Args
    ----
    filename : str
        Input filename.
    gene_id : str
        Input gene ID.

    """
    seqids, strands, locs = []
    with open(filename, "r") as input_handle:
        for line in input_handle:
            if line.startswith("#"):
                continue
            li = line.strip().split("\t")
            if len(li) != 9:
                continue
            seqid, src, kind, start, end, score, strand, phase, attributes = li

            if kind != "CDS":
                continue

            if gene_id in attributes:
                seqids.append(seqid)
                locs.append(int(start))
                locs.append(int(end))
                strands.append(strand)

    print(seqids[0], strands[0], min(locs), max(locs))