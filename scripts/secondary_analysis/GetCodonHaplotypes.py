
### Author: Andrew Valesano
### Project: Matlab Poliovirus Sequencing
### Purpose: Get haplotypes from BAM file.
### Working directory: Poliovirus_Intrahost

# Usage: python scripts/secondary_analysis/GetCodonHaplotypes.py --codon 2908 --outfile /Users/avalesano/Desktop/test.csv --infiles data/BAM_example/*.bam
# --codon is the first position of the codon to investigate

# Usage on GL:
# python GetCodonHaplotypes.py --codon 2908 --outfile codon_2908_frequencies.csv --infiles "plate*/data/aligned_output/removed_duplicates/*.bam"

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


def GetCodonFrequencies(codon_start_site, file):

    samfile = pysam.AlignmentFile(file, "rb")
    basename = file.split("/")[2]
    id = basename.split(".")[0]

    pos1 = codon_start_site + 0 - 1 # minus one because SAM positions are zero-based
    pos2 = codon_start_site + 1 - 1
    pos3 = codon_start_site + 2 - 1

    read_depth = 0
    reads = samfile.fetch("Sabin2ref", pos1, pos3)
    residue_counts = {}
    for read in reads:

        mapq_cutoff = 20
        if(read.mapq < mapq_cutoff):
            continue

        positions = read.get_reference_positions()
        if(pos1 not in positions or pos2 not in positions or pos3 not in positions):
            continue

        idx1 = positions.index(pos1)
        idx2 = positions.index(pos2)
        idx3 = positions.index(pos3)

        phred1 = ord(read.qual[idx1]) - 33
        phred2 = ord(read.qual[idx2]) - 33
        phred3 = ord(read.qual[idx3]) - 33

        phred_cutoff = 20
        if(phred1 < phred_cutoff or phred2 < phred_cutoff or phred3 < phred_cutoff):
            continue

        base1 = read.query_sequence[idx1]
        base2 = read.query_sequence[idx2]
        base3 = read.query_sequence[idx3]
        codon = Seq(base1 + base2 + base3)
        residue = str(codon.translate())
        read_depth = read_depth + 1
        residue_counts[residue] = residue_counts.get(residue, 0) + 1

    residue_counts_df = pd.DataFrame(list(residue_counts.items()), columns=['Residue', 'Count'])
    residue_counts_df = residue_counts_df.assign(Frequency = residue_counts_df.Count / read_depth)
    residue_counts_df = residue_counts_df.loc[residue_counts_df["Frequency"] > 0.01]
    residue_counts_df = residue_counts_df.assign(ID = id)
    residue_counts_df = residue_counts_df.assign(Depth = read_depth)

    return(residue_counts_df)


# ========================= Main =============================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--infiles', action="store", dest="infiles")
    parser.add_argument('--outfile', action = "store", dest = "outfile")
    parser.add_argument('--codon', action = "store", dest = "codon", type = int)
    args = parser.parse_args()

    all = pd.DataFrame()
    for file in glob.glob(args.infiles):
        filepath = Path(file)
        if(filepath.is_file() is False):
            continue
        results = GetCodonFrequencies(args.codon, file)
        all = pd.concat([all, results])

    print(all)
    #all.to_csv(args.outfile, index = False)


if __name__ == "__main__":
    main()
