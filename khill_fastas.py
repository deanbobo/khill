#!/usr/bin/env python3

import os
import sys
import argparse
import re
import mmh3
import math
import multiprocessing
from functools import partial
import logging


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
    parser.add_argument('-s', '--stats', help='File path for stats output', default=os.path.join(os.getcwd(), 'khill_out.txt'))
    parser.add_argument('-f', '--force', action='store_true', help='Force overwrite of the stats output file if it already exists.')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads to use. Default is 4.')
    parser.add_argument('-l', '--loglevel', help='Set the logging level (e.g., DEBUG, INFO, WARNING, ERROR, CRITICAL)', default='INFO')
    
    return parser.parse_args()


# logging setup
def setup_logging(args):
    logging.basicConfig(
        level=args.loglevel,  # Always log to console
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),  # Print to console
        ]
    )
    


# reverse compliment for canonical option
def revdnacomp(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))


# get the list of genomes for groups
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
                    logging.error(f"Error: Incorrect number of columns in line: {line.strip()}. Expected 2 columns.")
                    sys.exit(1)

                group_id, genome_path = parts
                if group_id not in group_genomes:
                    group_genomes[group_id] = []
                group_genomes[group_id].append(genome_path)

    # Process directory if provided
    if directory:
        # Use the directory name or path as the default group name
        # This can be customized based on your preference
        default_group_id = os.path.basename(directory) if os.path.basename(directory) else "default_group"
        group_genomes[default_group_id] = []
        for dir_path, _, file_names in os.walk(directory):
            for f in file_names:
                if re.match(r".*\.(fasta|fa|fna)$", f):
                    genome_path = os.path.join(dir_path, f)
                    group_genomes[default_group_id].append(genome_path)

    return group_genomes

# checking files
def check_input_files(groups):
    missing_files = []
    for group_id, genomes in groups.items():
        for genome_path in genomes:
            if not os.path.exists(genome_path):
                missing_files.append(genome_path)

    if missing_files:
        logging.error(f"{len(missing_files)} genome files are missing:")
        for file in missing_files:
            logging.error(f"Missing file: {file}")
        sys.exit(1)


# hash genomes
def process_single_genome(g, args, non_acgt_re):
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
                #the_hash = mmh3.hash64(cfrag, 460)[0]
                the_hash = mmh3.hash64(cfrag, 460, signed=False)[0]
                if the_hash < args.thresh:
                    HASHF.write(f'{the_hash}\n')
                kc += 1
                if kc >= args.max_kmers:
                    break
            i += args.motion


#  multiprocessing for hashing of genomes
def process_genomes(genomes, args):
    non_acgt_re = re.compile("[BDEFHIJKLMNOPQRSUVWXYZ0-9]", re.IGNORECASE)

    # Use multiprocessing to process genomes in parallel
    num_cores = min(len(genomes),args.threads)
    logging.info(f"Using {num_cores} cores for multiprocessing.")

    pool = multiprocessing.Pool(processes=num_cores)
    
    # Use functools.partial to pass additional arguments to process_single_genome
    process_func = partial(process_single_genome, args=args, non_acgt_re=non_acgt_re)

    # Process genomes in parallel
    pool.map(process_func, genomes)

    pool.close()
    pool.join()


# calculate kl
def kl(kl_g, kmers, kl_kpg, kl_sg, kl_allk):
    kl_val = 0
    for uniqmer in kmers.keys():

        # a particular kmer in a particular genome
        psi = 0
        if kl_g in kmers[uniqmer]:
            psi = kmers[uniqmer][kl_g] / kl_kpg[kl_g]

        # heuristic
        if psi == 0:
            continue

        # a particular kmer in the sampled pan genome
        pts = 0
        for subgs in kl_sg:
            if subgs in kmers[uniqmer]:
                pts += kmers[uniqmer][subgs]
        pss = pts / kl_allk
        
        # kl divergence
        kli = psi * (math.log(psi / pss))
        kl_val += kli

    return kl_val


# Function to calculate kl for a single genome and write betaent
def calculate_kl_for_genome(sg, kmers, kmersperg, genomes, allkmers, args):
    filebase = os.path.basename(sg)
    betaent_path = f"{args.output}/{filebase}.betaent"
    
    kl_val = kl(sg, kmers, kmersperg, genomes, allkmers) 
    
    # Calculate the weight (kmers per genome / all kmers)
    weight = kmersperg[sg] / allkmers
    
    be = kl_val * weight
    
    logging.info(f"Writing beta entropy for {filebase}")
    with open(betaent_path, 'w') as f:
        f.write(f"{filebase}\t{be}\t{kl_val}\t{weight}\n")
    
    return filebase 





def collate_hashes_write_betaent(genomes, args):
    # collate kmers
    kmers = {}
    allkmers = 0
    kmersperg = {}        
    for g in genomes:
        filebase = os.path.basename(g)
        logging.info(f"Collating hashes for {g}")
        try:
            with open(f"{args.output}/{filebase}.hashes", 'r') as HASHES:
                for h in HASHES:
                    h = h.strip()
                    h = int(h) #convert to integer
                    if h not in kmers:
                        kmers[h] = {}
                    if g not in kmers[h]:
                        kmers[h][g] = 0
                    # Increment the count for the k-mer in the genome
                    kmers[h][g] += 1
                    # Increment the total k-mers for the genome
                    if g not in kmersperg:
                        kmersperg[g] = 0
                    kmersperg[g] += 1
                    # Increment the total number of k-mers overall
                    allkmers += 1
        except FileNotFoundError:
            logging.warning(f"Warning: Hash file for {filebase} not found.")
            
    
    # Parallelize the calculation of kl and writing betaent
    num_threads = min(len(genomes), args.threads)
    logging.info(f"Using {num_threads} cores for parallel KL calculations.")
    
    # Set up a Pool for parallel processing
    pool = multiprocessing.Pool(processes=num_threads)
    
    # Partial function to pass the shared arguments (kmers, kmersperg, etc.)
    process_func = partial(calculate_kl_for_genome, kmers=kmers, kmersperg=kmersperg, genomes=genomes, allkmers=allkmers, args=args)
    
    # Run the calculations in parallel
    pool.map(process_func, genomes)
    
    pool.close()
    pool.join()


def calculate_khill(group_id, genomes, args):
    betaent = 0
    for sg in genomes:
        logging.info(f"Calculating KL for {sg}")
        filebase = os.path.basename(sg)
        try:
            with open(f"{args.output}/{filebase}.betaent", 'r') as HILL:
                for line in HILL:
                    parts = line.strip().split("\t")
                    # Check if the line contains 4 columns as expected
                    if len(parts) == 4:
                        # The second column (parts[1]) is the beta entropy (be)
                        betaent += float(parts[1])
                    else:
                        logging.warning(f"Unexpected format in {filebase}.betaent: {line}")
        except FileNotFoundError:
            logging.warning(f"Warning: Beta diversity file for {filebase} not found.")
    
    # Calculate beta diversity as the exponent of the accumulated betaent
    betadiv = math.exp(betaent)
    
    # Determine the stats file path
    stats_file_path = args.stats if args.stats else f"{args.output}/stats.txt"
    
    # Open stats file in append mode and write results
    with open(stats_file_path, 'a') as f:
        # Check if the file is empty to write the header
        if os.stat(stats_file_path).st_size == 0:
            f.write("#group_id\tbetaent\tbetadiv\n")
        
        # Write the group data (group ID, beta entropy, beta diversity)
        f.write(f"{group_id}\t{betaent}\t{betadiv}\n")


def process_groups(groups, args):
    for group_id, genomes in groups.items():
        logging.info(f"Processing group {group_id} with {len(genomes)} genomes.")
        
        # hash the genomes
        process_genomes(genomes, args)
        
        # collate hashes
        collate_hashes_write_betaent(genomes, args) 
        
        # from betaent files, calculate khill
        calculate_khill(group_id, genomes, args)


def main():
    args = parse_arguments()
    
    # Set up logging
    setup_logging(args)
    
    # Ensure the output directory exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)
        logging.info(f"Creating output directory: {args.output}")
    
    # Custom check for --directory or --file_table
    if not args.directory and not args.file_table:
        logging.error("Error: Either --directory or --file_table must be provided.")
        sys.exit(1)  # Exit the program with an error code
    
    # Check if the stats file already exists. Exit if it does, and --force not used
    stats_file_path = args.stats if args.stats else f"{args.output}/stats.txt"
    if os.path.isfile(stats_file_path) and not args.force:
        logging.error(f"Error: Stats file '{stats_file_path}' already exists. Use --force to overwrite.")
        sys.exit(1)
    elif os.path.isfile(stats_file_path) and args.force:
        logging.info(f"Overwriting stats file '{stats_file_path}'")
        os.remove(stats_file_path)
    
    # Calculate thresh based on hash_sketch
    maxhash = 2**64
    args.thresh = maxhash / args.hash_sketch
    
    # Get the sample groups
    groups = list_genomes(args.directory, args.file_table)
    
    # Check that all input files exist 
    check_input_files(groups)
    
    # Process each genome
    process_groups(groups, args)  # Generate hashes
    
    logging.info("K-Hill calculation completed.")

if __name__ == "__main__":
    main()
