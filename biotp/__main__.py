import argparse
import sys

from biotp import fasta
from biotp import gff
from biotp import jcvi
from biotp import json

def main():
    args = parse_args()
    function = functions[args.function]
    try:
        function(*args.args)
    except TypeError as e:
        print(e)
        sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('function', choices=functions)
    parser.add_argument('args', type=str, nargs='*')
    args = parser.parse_args()
    return args

functions = {
    'output_ids': fasta.output_ids,
    'rename_header': fasta.rename_header,
    'rename_headers_to_features': fasta.rename_headers_to_features,
    'slice_headers_by_ids': fasta.slice_headers_by_ids,
    'slice_seq_flanking_region': fasta.slice_seq_flanking_region,
    'slice_seq_upstream_region': fasta.slice_seq_upstream_region,
    'split_multi_into_single': fasta.split_multi_into_single,
    'merge_seqs_by_seqids': fasta.merge_seqs_by_seqids,
    'remove_n_from_seqs': fasta.remove_n_from_seqs,
    'slice_lines_by_seqids': gff.slice_lines_by_seqids,
    'output_seqid_strand_locs_by_proid': gff.output_seqid_strand_locs_by_proid,
    'extract_exon_intron': gff.extract_gene_exon_intron,
    'export_lengths_counts_of_gene_exon_intron': gff.export_lengths_counts_of_gene_exon_intron,
    'output_blocks': jcvi.output_blocks,
    'output_longest_one_to_two_microsynteny': jcvi.output_longest_one_to_two_microsynteny,
    'output_longest_one_to_three_microsynteny': jcvi.output_longest_one_to_three_microsynteny,
    'output_besthit_one_to_two_microsynteny': jcvi.output_besthit_one_to_two_microsynteny,
    'output_acc_org_asm': json.output_acc_org_asm
}

if __name__ == "__main__":
    main()



