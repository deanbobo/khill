# K-Hill

K-Hill is a method and software package that can quantify molecular diversity in pangenomes. The method draws on information theory (Shannon Diversity) to quantify richness and evenness of K-mers between genomes in groups of samples. The approach is computationally efficient - not relying on databases, alignments, or genome graphs. 

This `khill_fastas.py` script can be run two ways:
1. with --directory which  will measure K-Hill in a "pangeome" which includes all fasta files in the specified directory.
2. with --file_table which allows for grouping of genomes to be measured separately.

the file_table should have two columns. For example:

| #group_id | fasta_file_path |
| --- | --- |
| groupA | /path/to/file1.fa |
| groupA | /path/to/file2.fa |
| groupB | /path/to/file3.fa |
| groupB | /path/to/file4.fa |

The software will calculate K-Hill for each distinct group_id.

## USAGE

**-h, --help**      

>show the help message and exit

**-l, --loglevel

>Set the logging level (e.g., DEBUG, INFO, WARNING, ERROR, CRITICAL)', default='INFO'

**-t, --threads**
>number of threads to use. Default=4
  
**-d DIRECTORY, --directory DIRECTORY**
>Directory of genomes
>If --directory not used, --file_table must be provided.
  
**-o OUTPUT, --output OUTPUT**
 
>Directory for output of run. Default is current working directory. It is recommended to specify an empty directory
  
**-k KMER_SIZE, --kmer_size KMER_SIZE**
                 
>Kmer size
    
**-m MOTION, --motion MOTION**
                        
>Defined kmer overlap; 1 will move kmers along one nucleotide at a time
 
  
**-hs HASH_SKETCH, --hash_sketch HASH_SKETCH**
                        
>There are 2^64 hashes available for this hashing algorithm. 
    1 will look at every K-mer. 100 will examine 2^64/100 or 1 percent of K-mers. 
    100 is usually an acceptable value and saves computational time.
    
**-c, --cannon**
  
>Enable canonical kmers (reverse complement taken into account).

  
**-mk MAX_KMERS, --max_kmers MAX_KMERS**
  
>Limit the number of kmers to analyze
  
  
**-ft FILE_TABLE, --file_table FILE_TABLE**
>Table of files. Requires 2 columns. 
   First column is the the group_id, 
   second column is the path to the fasta file. 

>If --file_table is not provided, --directory must
                        be used.
 
**-s STATS, --stats STATS**
                       
>File path for stats output

  
**-f, --force**       
  
>Force overwrite of the stats output file if it already exists.

## Citation
https://genome.cshlp.org/content/34/10/1651.abstract


