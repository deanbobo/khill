#!/usr/bin/env python3

from os import walk
import sys
import getopt
import re
import string
import mmh3
from urllib.request import urlopen
import os

#from multiprocessing import Process, Manager, Lock
from math import log, exp

def usage():
    print('''
    USAGE: 
    -d <input_directory> : directory of genomes (optional, but either -d or -l must be provided)
    -o <output_directory> : directory with output of run
    -k <kmer_size> : kmer size
    -m <motion> : defined kmer overlap; m=1 will move kmers along on nuc at time and so on
    -n <number> : sampling rate
    -p <procs> : forks
    -c <cannon> : set to 1 if you want cannoical kmers (rev comp taken into account)
    -s <kmer_lim> : limit of kmers to analyze
    -l <filelist> : list of files (optional, but either -d or -l must be provided)
    ''')
    sys.exit()

def revdnacomp(st):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(st))


def log2(n):
    return log(n) / log(2)

def log10(n):
    return log(n) / log(10)

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



try:
    opts, args = getopt.getopt(sys.argv[1:], "d:o:k:p:m:n:c:s:l:")
except getopt.GetoptError as err:
    print(err)
    usage()

filelist = '/Users/dean/Dropbox/amnh/kmers/khill/filelist.txt'
#indir = '/Users/dean/Dropbox/amnh/kmers/20220630/'
outdir = '/Users/dean/Dropbox/amnh/kmers/khill/test_out/'
kmer = 19
motion = 1
number = 1
procs = 1
cannon = 1
kmerlim = 1e99

for o, a in opts:
    if o == "-d":
        indir = a
    elif o == "-o":
        outdir = a
    elif o == "-k":
        kmer = int(a)
    elif o == "-m":
        motion = int(a)
    elif o == "-n":
        number = a
    elif o == "-p":
        procs = a
    elif o == "-c":
        cannon = a
    elif o == "-s":
        kmerlim = int(a)
    elif o == "-l":
        filelist = a



try: indir
except: indir = None

try: filelist
except: filelist = None

#if (indir and filelist) is None or (indir or filelist) is not None:
#    raise ValueError('Can only specify -d <in_dir> OR -l <file_list>. not both.')





outdir=os.path.normpath(outdir)
if indir:
   indir=os.path.normpath(indir)

genomes = []



# list to store files name
if indir:
    res = []
    for (dir_path, dir_names, file_names) in walk(indir):
        res.extend(file_names)
    r = re.compile(".*\.fasta|.*\.fa|.*\.fna" )
    genomes = list(filter(r.match, res))

if filelist is not None:
    with open(filelist, 'r') as file:
        lines = file.readlines()
    genomes = [line.strip() for line in lines]

print(genomes)

maxhash = 2**64
thresh = maxhash/number

for g in genomes:
    filebase = os.path.basename(g) 
    HASHF=open(f"{outdir}/{filebase}.hashes", 'w')
    print("counting kmers for " + g)
    kc = 0
    signal = 0
    if indir:
        FASTA = open(f"{indir}/{g}", 'r')
    elif filelist:
        FASTA = open(f"{g}", 'r')

    for line in FASTA.readlines():
        line = line.strip()
        if line.startswith('>'):
            next
        else:
            seq = line.upper()
            cnter = 1
            while(len(seq) >= kmer):
                ffrag = seq[0:kmer]
                
                #disallow any string with unknown bases
                if re.search("[BDEFHIJKLMNOPQRSUVWXYZ0-9]", ffrag, re.IGNORECASE):
                    #print('Warning: Found non-ACGT character in fasta sequence.')
                    seq = seq[motion:] #NOT SURE ABOT THIS
                    cnter += motion
                else:
                    if cannon == 1:
                        rcfrag = revdnacomp(ffrag)
                        if ffrag < rcfrag:
                            cfrag = ffrag
                        else:
                            cfrag = rcfrag
                    else:
                        cfrag = ffrag
                    the_hash =  mmh3.hash64(cfrag, 460)[0]

                    if the_hash < thresh:
                        #print("we're here")
                        HASHF.write('%d\n' % the_hash)
                    kc=kc+1
                    seq = seq[motion:]
                    cnter += motion
                if kc >= kmerlim:
                    signal=signal+1
                if signal > 0:
                    break
        if signal > 0:
            break
    FASTA.close()
    HASHF.close()

# collate hashses
kmers = {}
allkmers = 0
kmersperg = {}
for g in genomes:
    filebase = os.path.basename(g)
    print(f"Collating hashes for {g}")
    HASHES = open(f"{outdir}/{filebase}.hashes", 'r')

    for h in HASHES.readlines():
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

# KLdiv
for sg in genomes:
    filebase = os.path.basename(sg)
    print(f"Calculating KL for {sg}")
    kl_val = kl(sg, kmers, kmersperg, genomes, allkmers)
    weight = (kmersperg[sg] / allkmers)
    be = kl_val * weight
    BETA = open(f"{outdir}/{filebase}.betaent", 'w')
    BETA.write(f"{sg}\t{be}\t{kl_val}\t{weight}\n")

# effective number of genomes (hill number for beta diversity)
betaent = 0
for sg in genomes:
    filebase = os.path.basename(sg)
    HILL=open(f"{outdir}/{filebase}.betaent", 'r')
    for line in HILL.readlines():
        line = line.strip()
        line = line.split("\t")
        betaent += float(line[1])
betadiv = exp(betaent)
with open(f"{outdir}/stats.txt", 'w') as f:
    f.write(f"betaEnt={betaent}\n")
    f.write(f"betaDiv={betadiv}\n")