import csv
import os

def slice_lines_by_seqids(input_filename, output_filename, *seqids):
    with open(input_filename, mode='r') as input_handle, open(output_filename, 'w') as output_handle:
        for line in input_handle:
            li = line.rstrip()
            if li.startswith('#'):
                continue  
            li = li.split('\t')
            seqid = li[0]
            if seqid not in seqids:
                continue
            output_handle.write(line)

def output_seqid_strand_locs_by_pepid(filename, pepid):
    seq_li = []
    strand_li = []
    loc_li = []
    with open(filename, mode='r') as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith('#'):
                continue  
            li = line.split('\t')
            items = li[8].split(';')
            for item in items:
                if item.startswith('protein_id='):
                    if pepid in item:
                        seq_li.append(str(li[0]))
                        loc_li.append(int(li[3]))
                        loc_li.append(int(li[4]))
                        strand_li.append(str(li[6]))
    seq = seq_li[0]
    strand = strand_li[0]
    start = min(loc_li)
    end = max(loc_li)
    print(seq, strand, start, end)

def generate_coordinate_all_introns(input_filename, output_filename):
    genes = {}
    
    with open(input_filename, 'r') as input_handle:
        for line in input_handle:
            if line.startswith('#'):
                continue
            li = line.strip().split('\t')
            if len(li) != 9:
                continue
            seqid, src, kind, start, end, score, strand, phase, attributes = li
            start, end = int(start), int(end)

            ## attributes hash
            attr_dict = {}
            for attr in attributes.split(';'):
                key, value = attr.split('=')
                attr_dict[key] = value

            ## gene line
            if kind == 'gene':
                gene_id = attr_dict['locus_tag']
                genes[gene_id] = {
                    'seq': seqid,
                    'src': src,
                    'start': start, 
                    'end': end, 
                    'strand': strand, 
                    'protein_id': set(),
                    'exons': []
                }

            ## exon line
            elif kind == 'exon':
                exon_id = attr_dict['locus_tag']
                if exon_id not in genes:
                    continue
                genes[exon_id]['exons'].append((start, end))

            ## CDS line
            elif kind == 'CDS':
                cds_id = attr_dict['locus_tag']
                if cds_id not in genes:
                    continue
                protein_id = attr_dict['protein_id']
                genes[cds_id]['protein_id'].add(protein_id)

    with open(output_filename, 'w') as output_handle:
        for gene_id, gene_info in genes.items():
            exons = sorted(gene_info['exons'])
            if len(gene_info['protein_id']) != 1:
                continue
            for protein_id in gene_info['protein_id']:
                intron_count = 1
                for j in range(len(exons) - 1):
                    intron_start = exons[j][1] + 1
                    intron_end = exons[j+1][0] - 1
                    if intron_start < intron_end:
                        intron_id = f"{protein_id}_intron_{intron_count}"
                        output_handle.write(f"{gene_info['seq']}\t{gene_info['src']}\tintron\t{intron_start}\t{intron_end}\t.\t{gene_info['strand']}\t.\tprotein_id={protein_id};intron_id={intron_id}\n")
                        intron_count += 1
