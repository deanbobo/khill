#!/usr/bin/env python3

from os import walk
import sys
import argparse
import re
import mmh3
import os
from math import log, exp


def parse_arguments():
    parser = argparse.ArgumentParser(description='Genome k-mer analysis tool.')
    parser.add_argument('-d', '--directory', help='Directory of genomes (optional, but either --directory or --file_table must be provided)', default=None)
    parser.add_argument('-o', '--output', help='Directory for output of run. Default is current working directory. It is recommended to specify an empty directory', default=os.getcwd())
    parser.add_argument('-k', '--kmer_size', type=int, default=19, help='Kmer size')
    parser.add_argument('-m', '--motion', type=int, default=1, help='Defined kmer overlap; 1 will move kmers along one nucleotide at a time')
    parser.add_argument('-hs', '--hash_sketch', type=int, default=1, help='There are 2^64 hashes available for this hashing algorithm. 1 will look at every K-mer. 100 will examine 2^64/100 or 1 percent of K-mers. 100 is usually an acceptable value and saves computational time.')
    parser.add_argument('-c', '--cannon', action='store_true', help='Enable canonical kmers (reverse complement taken into account).')
    parser.add_argument('-mk', '--max_kmers', type=float, default=1e99, help='Limit the number of kmers to analyze')
    parser.add_argument('-ft', '--file_table', help='Table of files. Requires 2 columns. First column is the the group_id, second column is the path to the fasta file. If --file_table is not provided, --directory must be used.')
    parser.add_argument('-s', '--stats', help='File path for stats output', default='/Users/dean/Dropbox/amnh/kmers/khill/test_out/stats.out')
    parser.add_argument('-f', '--force', action='store_true', help='Force overwrite of the stats output file if it already exists.')


    return parser.parse_args()

def revdnacomp(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))


def list_genomes(directory=None, file_table=None):
    group_genomes = {}

    # Process file_table if provided
    if file_table:
        with open(file_table, 'r') as file:
            for line in file:
                # Skip lines starting with '#'
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                # Check for the exact number of columns
                if len(parts) != 2:
                    print(f"Error: Incorrect number of columns in line: {line.strip()}. Expected 2 columns.")
                    sys.exit(1)

                group_id, genome_path = parts
                if group_id not in group_genomes:
                    group_genomes[group_id] = []
                group_genomes[group_id].append(genome_path)

    # Process directory if provided
    if directory:
        # Use the directory name or path as the default group name
        # This can be customized based on your preference
        default_group_id = path.basename(directory) if path.basename(directory) else "default_group"
        group_genomes[default_group_id] = []
        for dir_path, _, file_names in walk(directory):
            for f in file_names:
                if re.match(r".*\.(fasta|fa|fna)$", f):
                    genome_path = path.join(dir_path, f)
                    group_genomes[default_group_id].append(genome_path)

    return group_genomes

def check_input_files(groups):
    missing_files = []
    for group_id, genomes in groups.items():
        for genome_path in genomes:
            if not os.path.exists(genome_path):
                missing_files.append(genome_path)

    if missing_files:
        print("Error: The following genome files do not exist:")
        for file in missing_files:
            print(file)
        sys.exit(1)


def process_genomes(genomes, args):
    non_acgt_re = re.compile("[BDEFHIJKLMNOPQRSUVWXYZ0-9]", re.IGNORECASE)

    for g in genomes:
        filebase = os.path.basename(g)
        hash_file_path = f"{args.output}/{filebase}.hashes"
        
        with open(hash_file_path, 'w') as HASHF, open(g, 'r') as FASTA:
            seq = ""
            for line in FASTA:
                if not line.startswith('>'):
                    seq += line.strip().upper()
            
            kc = 0
            i = 0
            while i + args.kmer_size <= len(seq):
                ffrag = seq[i:i+args.kmer_size]
                if not non_acgt_re.search(ffrag):
                    cfrag = ffrag if not args.cannon else min(ffrag, revdnacomp(ffrag))
                    the_hash = mmh3.hash64(cfrag, 460)[0]
                    if the_hash < args.thresh:
                        HASHF.write(f'{the_hash}\n')
                    kc += 1
                    if kc >= args.max_kmers:
                        break
                i += args.motion


def kl(kl_g, kmers, kl_kpg, kl_sg, kl_allk):
    kl_val = 0
    for uniqmer in kmers.keys():

        # a particular kmer in a particular genome
        if kl_g not in kmers[uniqmer]:
            next
        else:
            psi = kmers[uniqmer][kl_g] / kl_kpg[kl_g]

            # heuristic
            if psi == 0:
                continue

            # a particular kmer in the sampled pan genome
            pts = 0
            for subgs in kl_sg:
                if subgs not in kmers[uniqmer]:
                    next
                else:
                    pts += kmers[uniqmer][subgs]
            pss = pts / kl_allk
            
            # kl divergence
            kli = psi * (log(psi/pss))
            kl_val += kli

    return kl_val

def collate_hashes_write_betaent(genomes, args):
    # collate kmers
    kmers = {}
    allkmers = 0
    kmersperg = {}        
    for g in genomes:
        filebase = os.path.basename(g)
        print(f"Collating hashes for {g}")
        try:
            with open(f"{args.output}/{filebase}.hashes", 'r') as HASHES:
                for h in HASHES:
                    h = h.strip()
                    h=int(h)
                    if h not in kmers:
                        kmers[h] = {}
                    if g not in kmers[h]:
                        kmers[h][g] = 0
                        kmers[h][g] += 1
                    if g not in kmersperg:
                        kmersperg[g] = 0
                    kmersperg[g] += 1
                    allkmers += 1
        except FileNotFoundError:
            print(f"Warning: Hash file for {filebase} not found.")

    # calculate kl and write betaent
    for sg in genomes:
        filebase = os.path.basename(sg)
        betaent_path = f"{args.output}/{filebase}.betaent"
    
        kl_val = kl(sg, kmers, kmersperg, genomes, allkmers) 
    
        with open(betaent_path, 'w') as f:
            f.write(f"{filebase}\t{kl_val}\n")


def calculate_khill(group_id, genomes, args):
    betaent = 0
    for sg in genomes:
        print(f"Calculating KL for {sg}")
        filebase = os.path.basename(sg)
        try:
            with open(f"{args.output}/{filebase}.betaent", 'r') as HILL:
                for line in HILL:
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        betaent += float(parts[1])
        except FileNotFoundError:
            print(f"Warning: Beta diversity file for {filebase} not found.")
    
    betadiv = exp(betaent)
    
    # Determine the stats file path outside of this function, if not already
    stats_file_path = args.stats if args.stats else f"{args.output}/stats.txt"
    
    # Open in 'a' mode to append if the file exists after handling --force in main
    with open(stats_file_path, 'a') as f:
        # Check if the file is empty to decide whether to write the header
        if os.stat(stats_file_path).st_size == 0:
            f.write("#group_id\tbetaent\tbetadiv\n")
        
        # Write the group data
        f.write(f"{group_id}\t{betaent}\t{betadiv}\n")


def process_groups(groups, args):
    for group_id, genomes in groups.items():
        print(f"Processing group {group_id} with {len(genomes)} genomes.")
        # Here you would adjust the functions called within this loop
        # to handle processing of genomes within the same group together.
        process_genomes(genomes, args)  # Adjust this function as necessary.
        # You might need to modify process_genomes and related functions
        # to handle genomes as a group, depending on how they're currently implemented.
        collate_hashes_write_betaent(genomes, args) # collate hashes
        calculate_khill(group_id, genomes, args)

def main():

    # get options
    args = parse_arguments()

    # Custom check for --directory or --file_table
    if not args.directory and not args.file_table:
        print("Error: Either --directory or --file_table must be provided.")
        sys.exit(1)  # Exit the program with an error code
    
    # check that the stats are going to print correctly
    stats_file_path = args.stats if args.stats else f"{args.output}/stats.txt"
    if os.path.isfile(stats_file_path) and not args.force:
        print(f"Error: Stats file '{stats_file_path}' already exists. Use --force to overwrite.")
        sys.exit(1)
    elif os.path.isfile(stats_file_path) and args.force:
        os.remove(stats_file_path)
    

    # Calculate thresh based on args.sampling_rate
    maxhash = 2**64
    args.thresh = maxhash / args.hash_sketch

    # get the sample groups
    groups = list_genomes(args.directory, args.file_table)

    # check that all input files exist 
    check_input_files(groups)

    # process each genome
    process_groups(groups, args) # generate hashes


if __name__ == "__main__":
    main()
