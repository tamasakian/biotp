#!/usr/bin/env python3

"""Library for processing GFF3 files.

Functions
---------
output_seqid_strand_locs_by_gene_id: Output seqid, strand, start, and end by partial match of gene id.
generate_upstream_regions: Generate upstream sequences for genes.
generate_downstream_regions: Generate downstream sequences for genes.

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
    seqids, strands, locs = [], [], []
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


def generate_upstream_regions(input_filename: str, output_filename: str, bp: int) -> None:
    """Generate upstream sequences for genes.

    Args
    ----
    input_filename : str
        Input GFF3 filename.
    output_filename : str
        Output GFF3 filename.
    bp : int
        Number of base pairs upstream to extract.

    """
    proteins = {}
    with open(input_filename, "r") as input_handle:
        for line in input_handle:
            if line.startswith("#"):
                continue
            li = line.strip().split("\t")
            if len(li) != 9:
                continue
            seqid, src, kind, start, end, score, strand, phase, attributes = li

            if kind != "CDS":
                continue
            start, end = int(start), int(end)

            attr_dict = {}
            for attr in attributes.split(";"):
                key, value = attr.split("=")
                attr_dict[key] = value
            if "protein_id" not in attr_dict:
                continue
            p = attr_dict["protein_id"]
            if p not in proteins:
                proteins[p] = {"seqid": seqid, "src": src, "strand": strand, "coordinates": []}
            proteins[p]["coordinates"].append((start, end))

    with open(output_filename, "w") as output_handle:
        bp = int(bp)
        for key, value in proteins.items():
            seqid = value["seqid"]
            src = value["src"]
            strand = value["strand"]
            start = min([coord[0] for coord in value["coordinates"]]) 
            end = max([coord[1] for coord in value["coordinates"]]) 

            if strand == "+":
                upstream_start = max(1, start - bp)
                upstream_end = start - 1
            else:
                upstream_start = end + 1
                upstream_end = end + bp

            output_handle.write(f"{seqid}\t{src}\tupstream\t{upstream_start}\t{upstream_end}\t.\t{strand}\t.\tprotein_id={key};upstream_id={key}_upstream_{bp}bp\n")


def generate_downstream_regions(input_filename: str, output_filename: str, bp: int) -> None:
    """Generate downstream sequences for genes.

    Args
    ----
    input_filename : str
        Input GFF3 filename.
    output_filename : str
        Output GFF3 filename.
    bp : int
        Number of base pairs downstream to extract.

    """
    proteins = {}
    with open(input_filename, "r") as input_handle:
        for line in input_handle:
            if line.startswith("#"):
                continue
            li = line.strip().split("\t")
            if len(li) != 9:
                continue
            seqid, src, kind, start, end, score, strand, phase, attributes = li

            if kind != "CDS":
                continue
            start, end = int(start), int(end)

            attr_dict = {}
            for attr in attributes.split(";"):
                key, value = attr.split("=")
                attr_dict[key] = value
            if "protein_id" not in attr_dict:
                continue
            p = attr_dict["protein_id"]
            if p not in proteins:
                proteins[p] = {"seqid": seqid, "src": src, "strand": strand, "coordinates": []}
            proteins[p]["coordinates"].append((start, end))

    with open(output_filename, "w") as output_handle:
        bp = int(bp)
        for key, value in proteins.items():
            seqid = value["seqid"]
            src = value["src"]
            strand = value["strand"]
            start = min([coord[0] for coord in value["coordinates"]]) 
            end = max([coord[1] for coord in value["coordinates"]]) 

            if strand == "+":
                downstream_start = end + 1
                downstream_end = end + bp
            else:
                downstream_start = max(1, start - bp)
                downstream_end = start - 1

            output_handle.write(f"{seqid}\t{src}\tdownstream\t{downstream_start}\t{downstream_end}\t.\t{strand}\t.\tprotein_id={key};downstream_id={key}_downstream_{bp}bp\n")
