
### Author: Andrew Valesano
### Project: Matlab Poliovirus Sequencing
### Purpose: For the donor population in pair 702, get as much info as we can on linkage of minor variants.
### Working directory: Poliovirus_Intrahost

# Usage: python scripts/secondary_analysis/GetCodonHaplotypes.py --codon 2908 --outfile /Users/avalesano/Desktop/test.csv --infiles data/BAM_example/*.bam

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

def GetReadCounts(pos1, pos2, samfile, var1, var2):

    pos1_samfile = pos1 - 1 # minus one because SAM positions are zero-based
    pos2_samfile = pos2 - 1

    if(pos1 > pos2):
        reads = samfile.fetch("Sabin2ref", pos2_samfile, pos1_samfile)
    elif(pos1 == pos2):
        reads = samfile.fetch("Sabin2ref", pos1_samfile, pos1_samfile + 1)
    else:
        reads = samfile.fetch("Sabin2ref", pos1_samfile, pos2_samfile)

    total_reads = 0
    linked_variant_reads = 0
    var1_reads = 0
    var2_reads = 0
    for read in reads:

        # Filter if both positions aren't in the read
        positions = read.get_reference_positions()
        if(pos1_samfile not in positions or pos2_samfile not in positions):
            continue

        # Filter if the MapQ is bad
        mapq_cutoff = 20
        if(read.mapq < mapq_cutoff):
            continue

        idx1 = positions.index(pos1_samfile)
        idx2 = positions.index(pos2_samfile)

        # Filter if the Phred scores are bad
        phred1 = ord(read.qual[idx1]) - 33
        phred2 = ord(read.qual[idx2]) - 33
        phred_cutoff = 20
        if(phred1 < phred_cutoff or phred2 < phred_cutoff):
            continue

        # Ok, count this read and record data
        total_reads = total_reads + 1

        base1 = read.query_sequence[idx1]
        base2 = read.query_sequence[idx2]

        if(base1 == var1):
            var1_reads = var1_reads + 1
        if(base2 == var2):
            var2_reads = var2_reads + 1
        if(base1 == var1 and base2 == var2):
            linked_variant_reads = linked_variant_reads + 1

    if(total_reads == 0):
        print(str(pos1) + "," + str(pos2) + ",0,0,0,0")
    else:
        freq1 = round(float(var1_reads) / total_reads, 4)
        freq2 = round(float(var2_reads) / total_reads, 4)
        freq_linked = round(float(linked_variant_reads) / total_reads, 4)
        print(str(pos1) + "," + str(pos2) + "," + str(total_reads) + "," + str(freq1) + "," + str(freq2) + "," + str(freq_linked))

    return()

# ========================= Main =============================

def main():

    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    polymorphic_positions = pd.read_csv("data/processed/variants_donor_polymorphic_702_minor.csv", index_col = None, header = 0, dtype = object)
    samfile = pysam.AlignmentFile("data/raw/plate3/131797.removed.bam", "rb")

    print("Position1,Position2,TotalReads,Freq1,Freq2,FreqLinked")
    for pos1 in polymorphic_positions["pos"]:
        for pos2 in polymorphic_positions["pos"]:

            # Get expected variant nucleotide at each site
            polymorphic_positions_1 = polymorphic_positions.loc[polymorphic_positions["pos"] == pos1]
            polymorphic_positions_2 = polymorphic_positions.loc[polymorphic_positions["pos"] == pos2]
            var1 = polymorphic_positions_1.iloc[0]['var']
            var2 = polymorphic_positions_2.iloc[0]['var']
            GetReadCounts(int(pos1), int(pos2), samfile, var1, var2)


if __name__ == "__main__":
    main()
