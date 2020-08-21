
### Author: Andrew Valesano
### Project: Matlab Poliovirus Sequencing
### Purpose: Get haplotypes from BAM file.
### Working directory: Poliovirus_Intrahost

# Usage: python scripts/secondary_analysis/GetCodonHaplotypes.py --codon 2908 --outfile /Users/avalesano/Desktop/test.csv --infiles data/BAM_example/*.bam
# --codon is the first position of the codon to investigate

# Usage on GL:
# python GetBaseCounts.py --base 481 --outfile base_481_frequencies.csv --infiles "plate*/data/aligned_output/removed_duplicates/*.bam"

# ======================= Import modules ======================

import argparse
import glob
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from pathlib2 import Path
import pysam
from Bio.Seq import Seq

# ========================= Functions =====================


def GetBaseCounts(site, file):

    samfile = pysam.AlignmentFile(file, "rb")
    basename = file.split("/")[4]
    id = basename.split(".")[0]
    print(id)

    pos = site - 1 # minus one because SAM positions are zero-based

    read_depth = 0
    reads = samfile.fetch("Sabin2ref", pos)
    base_counts = {}
    for read in reads:

        mapq_cutoff = 20
        if(read.mapq < mapq_cutoff):
            continue

        positions = read.get_reference_positions()
        if(pos not in positions):
            continue

        idx = positions.index(pos)

        phred = ord(read.qual[idx]) - 33

        phred_cutoff = 20
        if(phred < phred_cutoff):
            continue

        base = read.query_sequence[idx]
        read_depth = read_depth + 1
        base_counts[base] = base_counts.get(base, 0) + 1

    base_counts = pd.DataFrame(list(base_counts.items()), columns=['Base', 'Count'])
    base_counts = base_counts.assign(Frequency = base_counts.Count / read_depth)
    base_counts = base_counts.loc[base_counts["Frequency"] > 0.005]
    base_counts = base_counts.assign(ID = id)
    base_counts = base_counts.assign(Depth = read_depth)

    return(base_counts)


# ========================= Main =============================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--infiles', action="store", dest="infiles")
    parser.add_argument('--outfile', action = "store", dest = "outfile")
    parser.add_argument('--base', action = "store", dest = "base", type = int)
    args = parser.parse_args()

    all = pd.DataFrame()
    for file in glob.glob(args.infiles):
        filepath = Path(file)
        if(filepath.is_file() is False):
            continue
        print("Working on file: " + file)
        results = GetBaseCounts(args.base, file)
        all = pd.concat([all, results])

    print(all)
    all.to_csv(args.outfile, index = False)


if __name__ == "__main__":
    main()
