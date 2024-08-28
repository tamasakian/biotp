import argparse
import sys

from biotp import blast
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
    'slice_rgo_by_hgt': blast.slice_rgo_by_hgt,
    'detect_hgt': blast.detect_hgt, 
    'slice_results_by_rgo': blast.slice_results_by_rgo,
    'slice_results_by_rgo_2': blast.slice_results_by_rgo_2, 
    'slice_hits_by_crossover_group': blast.slice_hits_by_crossover_group,
    'output_besthit_for_subgroup': blast.output_besthit_for_subgroup,
    'output_ids': fasta.output_ids,
    'rename_header': fasta.rename_header,
    'rename_headers_to_features': fasta.rename_headers_to_features,
    'rename_headers_feature': fasta.rename_headers_feature,
    'prefix_to_headers': fasta.prefix_to_headers,
    'slice_records_by_names': fasta.slice_records_by_names, 
    'slice_headers_by_ids': fasta.slice_headers_by_ids,
    'slice_seq_flanking_region': fasta.slice_seq_flanking_region,
    'slice_seq_upstream_region': fasta.slice_seq_upstream_region,
    'split_multi_into_single': fasta.split_multi_into_single,
    'merge_seqs_by_seqids': fasta.merge_seqs_by_seqids,
    'remove_n_from_seqs': fasta.remove_n_from_seqs,
    'generate_all_introns': fasta.generate_all_introns,
    'slice_lines_by_seqids': gff.slice_lines_by_seqids,
    'output_seqid_strand_locs_by_pepid': gff.output_seqid_strand_locs_by_pepid,
    'generate_coordinate_all_introns': gff.generate_coordinate_all_introns,
    'make_dict_pepid': gff.make_dict_pepid,
    'output_blocks': jcvi.output_blocks,
    'output_longest_one_to_one_microsynteny': jcvi.output_longest_one_to_one_microsynteny,
    'output_longest_one_to_two_microsynteny': jcvi.output_longest_one_to_two_microsynteny,
    'output_longest_one_to_three_microsynteny': jcvi.output_longest_one_to_three_microsynteny,
    'output_besthit_one_to_two_microsynteny': jcvi.output_besthit_one_to_two_microsynteny,
    "output_besthit_one_to_two_synteny": jcvi.output_besthit_one_to_two_synteny, 
    'output_acc_org_asm': json.output_acc_org_asm
}

if __name__ == "__main__":
    main()



