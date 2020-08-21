
### Author: Andrew Valesano
### Project: Matlab Poliovirus Sequencing
### Purpose: Get nucleotides at specific positions from a cVDPV alignment.
### Working directory: Poliovirus_Intrahost

# Usage for Stern et al. alignments:
# Usage for 5UTR: python ./scripts/secondary_analysis/GetBaseFromVDPV.py --infile data/reference/VDPV_alignments/Stern_5UTR_alignment_muscle.fasta --variants data/reference/VDPV_alignments/5UTR_variants.csv
# Usage for P1: python ./scripts/secondary_analysis/GetBaseFromVDPV.py --infile data/reference/VDPV_alignments/Stern_P1_alignment_sed_muscle.fasta --variants data/reference/VDPV_alignments/P1_variants.csv

# Usage for Famulare/IDM alignments:
# 5_UTR: python ./scripts/secondary_analysis/GetBaseFromVDPV.py --infile data/reference/IDM_Sabin2_Alignments/Sabin2.nt5NCR.alignment.final.fasta --variants data/reference/VDPV_alignments/5UTR_variants.csv
# P1: python ./scripts/secondary_analysis/GetBaseFromVDPV.py --infile data/reference/IDM_Sabin2_Alignments/Sabin2.ntP1.alignment.final.fasta --variants data/reference/VDPV_alignments/P1_variants.csv
# P2: python ./scripts/secondary_analysis/GetBaseFromVDPV.py --infile data/reference/IDM_Sabin2_Alignments/Sabin2.ntP2.alignment.final.fasta --variants data/reference/VDPV_alignments/P2_variants.csv
# P3: python ./scripts/secondary_analysis/GetBaseFromVDPV.py --infile data/reference/IDM_Sabin2_Alignments/Sabin2.ntP3.alignment.final.fasta --variants data/reference/VDPV_alignments/P3_variants.csv

# ======================= Import modules ======================

import argparse
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from pathlib2 import Path

# ========================= Functions =====================

def GetAlignmentProportions(pos, var, file, accessions):

    # Changed for IDM alignments
    if("5UTR" in file or "5NCR" in file): # offset = 1221 for P1 file; 0 for 5UTR
        offset = 0
    elif("P1" in file):
        offset = 747
    elif("P2" in file):
        offset = 3384
    elif("P3" in file):
        offset = 5109
    elif("3NCR" in file):
        offset = 7368
    else:
        print("Woah, some serious error occurred!")
        return

    position = pos - offset

    total = 0
    num_var = 0
    genomes_queried = 0
    for record in SeqIO.parse(file, "fasta"):
        id = record.id

        if(id not in accessions):
            continue
        seq = record.seq
        base = seq[position - 1] # starts at base 0

        if("Sabin2" in id):
            print("found it, but we shouldn't have!")

        if(base is "-"):
            continue

        total = total + 1
        if(base is var):
            num_var = num_var + 1

    print(str(total))
    return(round(float(num_var) / total, 2))


# ========================= Main =============================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', action="store", dest="infile")
    parser.add_argument('--variants', action="store", dest="variantfile")
    args = parser.parse_args()

    variants = pd.read_csv(args.variantfile, index_col = None, header = 0, dtype = object)

    with open("data/reference/IDM_Sabin2_Alignments/accessionsToUse.txt") as f:
        accessions = f.read().splitlines()

    print("Reference,Position,Variant,Proportion,Count,Type,Location")
    for index, row in variants.iterrows():
        proportion = GetAlignmentProportions(int(row["Position"]), row["Variant"], args.infile, accessions)
        print(row["Reference"] + "," + row["Position"] + "," + row["Variant"] + "," + str(proportion) + "," + row["Count"] + "," + row["Type"] + "," + row["Location"])


if __name__ == "__main__":
    main()
