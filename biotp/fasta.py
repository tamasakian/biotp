import re
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def output_ids(filename):
    with open(filename, 'r') as handle:
        for i, record in enumerate(SeqIO.parse(handle, 'fasta'), 0):
            id = record.id
            print(i, id)

def rename_header(input_filename, output_filename, new_id, new_name, new_description):
    with open(input_filename, mode='r') as input_handle:
        record = SeqIO.read(input_handle, 'fasta')
        ## Rewrite
        record.id = new_id
        record.name = new_name
        record.description = new_description
    with open(output_filename, mode='w') as output_handle:
        SeqIO.write(record, output_handle, 'fasta')

def rename_headers_to_features(input_filename, output_filename, feature):
    with open(input_filename, 'r') as input_handle, open(output_filename, 'w') as output_handle:
        for record in SeqIO.parse(input_handle, 'fasta'):
            match = re.search(r'\[' + feature + r'=([^\]]+)]', record.description)  
            if not match:
                match = re.search(r'\[' + feature + r'=([^\]]+)]', record.name)
                if not match:
                    match = re.search(r'\[' + feature + r'=([^\]]+)]', record.id)
                    if not match:
                        SeqIO.write(record, output_handle, 'fasta')
                        continue
            new_id = match.group(1)
            ## Rewrite
            record.id = new_id
            record.name = ""
            record.description = ""
            SeqIO.write(record, output_handle, 'fasta')

def slice_headers_by_ids(input_filename, output_filename, *ids):
    with open(input_filename, 'r') as input_handle, open(output_filename, 'w') as output_handle:
        for record in SeqIO.parse(input_handle, 'fasta'):
            id = record.id.split("lcl|")[1].split("_cds_")[0]
            if id not in ids:
                continue
            SeqIO.write(record, output_handle, 'fasta')

def split_multi_into_single(input_filename, output_dirname):
    with open(input_filename, 'r') as input_handle:
        for record in SeqIO.parse(input_handle, 'fasta'):
            output_id = record.id
            output_filename = f'{output_dirname}/{output_id}.fasta'
            SeqIO.write(record, output_filename, 'fasta')

def merge_seqs_by_seqids(input_filename, output_filename):
    merged_seqs = defaultdict(str)
    with open(input_filename, 'r') as input_handle:
        for record in SeqIO.parse(input_handle, 'fasta'):
            merged_seqs[record.id] += str(record.seq)
    with open(output_filename, 'w') as output_handle:
        seq_records = []
        for seqid, sequence in merged_seqs.items():
            seq_record = SeqRecord(Seq(sequence), id=seqid, description='')
            seq_records.append(seq_record)
        SeqIO.write(seq_records, output_handle, 'fasta')

def slice_seq_flanking_region(input_filename, output_filename, output_id, output_name, output_description, strand, start, end, bp):
    with open(input_filename, mode='r') as input_handle:
        start = int(start)
        end = int(end)
        bp = int(bp)
        record = SeqIO.read(input_handle, 'fasta')
        if strand == '+':
            sliced_seq = record.seq[start-1-bp:end-1+bp]
        elif strand == '-':
            sliced_seq = record.seq[start-3-bp:end+bp].reverse_complement()
        else:
            sliced_seq = None
        scodon = sliced_seq[bp:bp+3]
        length = len(sliced_seq)
        print(f'start codon: {scodon}, length: {length}')
    with open(output_filename, mode='w') as output_handle:
        record = SeqRecord(
            seq         = sliced_seq,
            id          = output_id,
            name        = output_name,
            description = output_description)
        SeqIO.write(record, output_handle, 'fasta')

def slice_seq_upstream_region(input_filename, output_filename, output_id, output_name, output_description, strand, start, end, bp):
    with open(input_filename, mode='r') as input_handle:
        start = int(start)
        end = int(end)
        bp = int(bp)
        record = SeqIO.read(input_handle, 'fasta')
        if strand == '+':
            sliced_seq = record.seq[start-1-bp:start+2]
        elif strand == '-':
            sliced_seq = record.seq[end-3:end+bp].reverse_complement()
        else:
            sliced_seq = None
        scodon = sliced_seq[bp:bp+3]
        length = len(sliced_seq)
        print(f'start-codon: {scodon}, length: {length}')
    with open(output_filename, mode='w') as output_handle:
        record = SeqRecord(
            seq         = sliced_seq,
            id          = output_id,
            name        = output_name,
            description = output_description)
        SeqIO.write(record, output_handle, 'fasta')

def remove_n_from_seqs(input_filename, output_filename):
    output_records = []
    with open(input_filename, mode='r') as input_handle:
        for record in SeqIO.parse(input_handle, 'fasta'):
            removed_seq = str(record.seq).replace('N', '')
            output_record = SeqRecord(
                seq = Seq(removed_seq),
                id = record.id,
                description = record.description)
            output_records.append(output_record)
    with open(output_filename, mode='w') as output_handle:
        SeqIO.write(output_records, output_handle, 'fasta')


