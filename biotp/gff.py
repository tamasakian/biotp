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

def output_seqid_strand_locs_by_proid(filename, proid):
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
                    if proid in item:
                        seq_li.append(str(li[0]))
                        loc_li.append(int(li[3]))
                        loc_li.append(int(li[4]))
                        strand_li.append(str(li[6]))
    seq = seq_li[0]
    strand = strand_li[0]
    start = min(loc_li)
    end = max(loc_li)
    print(seq, strand, start, end)

## 建設中
def extract_gene_exon_intron(input_filename, output_filename):
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
            attr_dict = {}
            for attr in attributes.split(';'):
                key, value = attr.split('=')
                attr_dict[key] = value
            if kind == 'gene':
                gene_id = attr_dict['ID']
                genes[gene_id] = {
                    'seq': seqid,
                    'src': src,
                    'start': start, 
                    'end': end, 
                    'strand': strand, 
                    'mRNAs': [],
                    'exons': []
                }
            elif kind == 'mRNA':
                parent_gene = attr_dict['Parent']
                genes[parent_gene]['mRNAs'].append((start, end, attr_dict['ID']))
            elif kind == 'exon':
                parent_mRNA = attr_dict['Parent']
                for gene_id, gene_info in genes.items():
                    for mRNA in gene_info['mRNAs']:
                        if mRNA[2] == parent_mRNA:
                            genes[gene_id]['exons'].append((start, end, parent_mRNA))
                            break
    with open(output_filename, 'w') as output_handle:
        for gene_id, gene_info in genes.items():
            ## gene
            output_handle.write(f'{gene_info['seq']}\t{gene_info['src']}\tgene\t{gene_info['start']}\t{gene_info['end']}\t.\t{gene_info['strand']}\t.\tID={gene_id}\n')
            mRNAs = gene_info['mRNAs']
            for mRNA in mRNAs:
                mRNA_start, mRNA_end, mRNA_id = mRNA
                ## mRNA
                output_handle.write(f'{gene_info['seq']}\t{gene_info['src']}\tmRNA\t{mRNA_start}\t{mRNA_end}\t.\t{gene_info['strand']}\t.\tID={mRNA_id};Parent={gene_id}\n')
                exons = [exon for exon in gene_info['exons'] if exon[2] == mRNA_id]
                exons.sort()
                for i, exon in enumerate(exons):
                    exon_start, exon_end, exon_parent = exon
                    ## exon
                    output_handle.write(f"{gene_info['seq']}\t{gene_info['src']}\texon\t{exon_start}\t{exon_end}\t.\t{gene_info['strand']}\t.\tID={exon_parent}_exon{i+1};Parent={exon_parent}\n")
                for j in range(len(exons) - 1):
                    intron_start = exons[j][1] + 1
                    intron_end = exons[j+1][0] - 1
                    if intron_start < intron_end:
                        ## intron
                        output_handle.write(f"{gene_info['seq']}\t{gene_info['src']}\tintron\t{intron_start}\t{intron_end}\t.\t{gene_info['strand']}\t.\tID={mRNA_id}_intron{j+1};Parent={mRNA_id}\n")

## 建設中
def export_lengths_counts_of_gene_exon_intron(input_filename, output_filename):
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
            attr_dict = {}
            for attr in attributes.split(';'):
                key, value = attr.split('=')
                attr_dict[key] = value
            if kind == 'gene':
                gene_id = attr_dict['ID']
                genes[gene_id] = {'exons': [], 'introns': [], 'mRNAs': []}
            elif kind == 'mRNA':
                parent_gene = attr_dict['Parent']
                genes[parent_gene]['mRNAs'].append((start, end, attr_dict['ID']))
            elif kind == 'exon':
                parent_mRNA = attr_dict['Parent']
                for gene_id, gene_info in genes.items():
                    if any(mRNA for mRNA in gene_info.get('mRNAs', []) if mRNA[2] == parent_mRNA):
                        genes[gene_id]['exons'].append((start, end))
                        break
            elif kind == 'intron':
                parent_mRNA = attr_dict['Parent']
                for gene_id, gene_info in genes.items():
                    if any(mRNA for mRNA in gene_info.get('mRNAs', []) if mRNA[2] == parent_mRNA):
                        genes[gene_id]['introns'].append((start, end))
                        break
    ## calculate lengths and counts
    result = []
    for gene_id, gene_info in genes.items():
        exon_length = sum(end - start + 1 for start, end in gene_info['exons'])
        intron_length = sum(end - start + 1 for start, end in gene_info['introns'])
        exon_count = len(gene_info['exons'])
        intron_count = len(gene_info['introns'])
        result.append({'gene': gene_id, 'kind': 'exon', 'length': exon_length, 'count': exon_count})
        result.append({'gene': gene_id, 'kind': 'intron', 'length': intron_length, 'count': intron_count})
    ## write to csv
    write_header = not os.path.exists(output_filename)
    with open(output_filename, 'a', newline='') as output_handle:
        fieldnames = ['gene', 'kind', 'length', 'count']
        writer = csv.DictWriter(output_handle, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        for row in result:
            writer.writerow(row)
