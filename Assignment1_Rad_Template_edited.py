#!/usr/bin/env python3
import argparse,pdb
import pyfaidx
from Bio.Restriction.Restriction_Dictionary import rest_dict
from typing import List, Tuple


# Task 1: Parsing Command-Line Arguments

def my_parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='This script performs either SingleRad or ddRad sequencing on a DNA FASTA file and analyzes variation from a VCF file.'
    )
    parser.add_argument('GenomeFile', type=str, help='Path to the genome FASTA file')
    parser.add_argument('VCFFile', type=str, help='Path to the VCF file')
    parser.add_argument('Mode', type=str, choices=['SingleRad', 'ddRad'], help='Sequencing mode (SingleRad or ddRad)')
    parser.add_argument('RE1', type=str, help='First restriction enzyme name')
    parser.add_argument('--RE2', type=str, required=False, help='Second restriction enzyme name (required for ddRad)')
    parsed_args = parser.parse_args()
    return parsed_args


# Task 2: Reading the FASTA File

def read_fasta(path_to_fasta: str, chromosome: str = 'NC_036780.1') -> str:
    fasta = pyfaidx.Fasta(path_to_fasta)
    dna = str(fasta[chromosome][:])
    return dna


# Task 3: Finding Motifs

def find_motifs(dna: str, motif: str) -> List[int]:
    positions = []
    start = 0
    while True:
        index = dna.find(motif, start)
        if index == -1:
            break
        positions.append(index)
        start = index + len(motif)
    return positions




# Task 4: Implementing SingleRad
def run_single_rad(dna: str, re1: str) -> List[Tuple[int, int]]:
    seq_length = 100  
    if re1 not in rest_dict:
        raise ValueError(f"Restriction enzyme {re1} not found in rest_dict.")
    motif = rest_dict[re1]['site']
    cut_sites = find_motifs(dna, motif)
    sequenced_sites = [(site - seq_length, site + seq_length) for site in cut_sites]
    return sequenced_sites




# Task 5: Implementing ddRad
def run_ddrad(dna: str, re1: str, re2: str) -> List[Tuple[int, int]]:
    min_size = 300
    max_size = 700
    seq_length = 100  # 100 bp fragment length for ddRad
    if re1 not in rest_dict or re2 not in rest_dict:
        raise ValueError("Restriction enzyme not found in rest_dict.")
    motif1 = rest_dict[re1]['site']
    motif2 = rest_dict[re2]['site']
    re1_sites = find_motifs(dna, motif1)
    re2_sites = find_motifs(dna, motif2)
    dd_sites = []
  
    for r1 in re1_sites:
        left_hits = [r2 for r2 in re2_sites if r2 < r1]
        left_hits_r1 = [x for x in re1_sites if x < r1]
        if left_hits:
            nearest_r2 = left_hits[-1] 
            if min_size < (r1 - nearest_r2) < max_size:
                if nearest_r2 > left_hits_r1[-1]:
                    dd_sites.append((nearest_r2, nearest_r2 + seq_length))  
                    dd_sites.append((r1 - seq_length, r1))

    for r1 in re1_sites:
        right_hits = [r2 for r2 in re2_sites if r2 > r1]
        right_hits_r1 = [x for x in re1_sites if x > r1]
        if right_hits:
            nearest_r2 = right_hits[0] 
            if min_size < (nearest_r2 - r1) < max_size:
                if nearest_r2 < right_hits_r1[0]:
                    dd_sites.append((r1, r1 + seq_length))  
                    dd_sites.append((nearest_r2 - seq_length, nearest_r2))  

    return dd_sites




# Task 6: Finding Variable Sites

import vcfpy

def find_variable_sites(vcf_file_path: str, sequenced_sites: list[tuple[int, int]]) -> list[tuple[int, int]]:
    with open(vcf_file_path, "r") as f:
        reader = vcfpy.Reader(f)
        sequenced_sites_variable = []
        for record in reader:
            if record.CHROM != 'NC_036780.1':
                break
            if len(record.calls) != 2:
                continue
            gt1 = record.calls[0].data.get('GT')
            gt2 = record.calls[1].data.get('GT')
            if (gt1 == '0/0' and gt2 == '1/1') or (gt1 == '1/1' and gt2 == '0/0'):
                matches = [site for site in sequenced_sites if site[0] < record.POS < site[1]]
                if matches:
                    sequenced_sites_variable.extend(matches)
        return sequenced_sites_variable
    
### MAIN CODE ###
# Once you think your code is complete, comment out the above testing section,
# and run the script.

def main():
    args = my_parse_args()
    target_chromosome = 'NC_036780.1'
    dna_string = read_fasta(path_to_fasta=args.GenomeFile, chromosome=target_chromosome)
    
args = my_parse_args()
target_chromosome = 'NC_036780.1'
dna_string = read_fasta(path_to_fasta=args.GenomeFile, chromosome=target_chromosome)

if args.Mode == 'SingleRad':
    print(f'running SingleRad for chromosome {target_chromosome} and restriction enzyme {args.RE1}')
    seq_sites = run_single_rad(dna_string, args.RE1)
elif args.Mode == 'ddRad':
    print(f'running ddRad for chromosome {target_chromosome} and restriction enzymes {args.RE1} and {args.RE2}')
    seq_sites = run_ddrad(dna_string, args.RE1, args.RE2)

input_len = len(dna_string)
n_sites = len(seq_sites)
sequencing_length = seq_sites[0][1] - seq_sites[0][0]
sequenced_fraction = sequencing_length * n_sites / input_len
print(f'{args.Mode} sequencing complete.')
print(f'\tlength of input sequence: {input_len}')
print(f'\tnumber of sequencing sites located: {n_sites}')
print(f'\tpercentage of nucleotides sequenced: {100 * sequenced_fraction:.3f}%')

print('checking for variation within sequenced sites')
variable_sites = find_variable_sites(args.VCFFile, sequenced_sites=seq_sites)
n_var_sites = len(variable_sites)
var_fraction = n_var_sites / n_sites
print(f'\tnumber of sites with variation: {n_var_sites}')
print(f'\tfraction of sites with variation: {var_fraction}')

print('\nAnalysis Complete')
