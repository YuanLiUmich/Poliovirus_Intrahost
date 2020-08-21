

### Author: Andrew Valesano
### Project: Matlab Poliovirus Sequencing
### Purpose: For samples with duplicates, use the consensus info from both replicates to get a complete consensus sequence.
### Use this after doing FilterConsensus.py.
### Working directory: Poliovirus_Intrahost

### Usage: python ./scripts/secondary_analysis/ResolveConsensusDuplicates.py --indir data/processed/consensus_filtered_final/ --outdir data/processed/consensus_resolved_final/ --meta data/metadata/specimens_Sabin2Positive_CTbelow37_processedFinal_deID_runInfo.csv --coverage data/raw_final/ > data/processed/consensus_resolved_final/conflictOutput.csv

# ======================= Import modules ======================

import argparse
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from pathlib2 import Path

# ========================= Functions =============================

def LoadCoverage(indir):

    coverage_files = "%s*/all.coverage.csv" % indir
    cov = []
    for file in glob.glob(coverage_files):
        print(file)
        c = pd.read_csv(file, index_col = None, header = 0, dtype = object)
        cov.append(c)
    coverage = pd.concat(cov, axis=0, ignore_index=True)
    coverage.columns = coverage.columns.astype(str)

    return coverage

def ResolveFasta(seqrecord1, seqrecord2, coverage, id1, id2, bio_id):

    for entry in seqrecord1:
        seq1 = entry.seq
    for entry in seqrecord2:
        seq2 = entry.seq

    genomeLength1 = len(seq1)
    genomeLength2 = len(seq2)
    if(genomeLength1 != genomeLength2):
        print("Woah the sequences should be the same length. Abort!")
        return None

    new_sequence = list("x" * genomeLength1)
    for i in range(0, genomeLength1, 1):
        base1 = seq1[i]
        base2 = seq2[i]
        new_base = "x"
        if(base1 is "N" and base2 is "N"):
            new_base = "N"
        elif(base1 is not "N" and base2 is "N"):
            new_base = base1
        elif(base1 is "N" and base2 is not "N"):
            new_base = base2
        elif(base1 is not "N" and base2 is not "N"):
            new_base = base1
            if(base1 is not base2):
                #new_base = "N" # option for the leaveN version
                print(bio_id + "," + str(i+1))
                #print("base1: " + base1 + "; base2: " + base2)
                #sample1_coverage = coverage.loc[coverage["Id"] == id1]
                #sample1_coverage = sample1_coverage.loc[sample1_coverage["chr"] == "Sabin2ref"]
                #sample1_coverage = sample1_coverage.loc[sample1_coverage["chr.pos"] == str(i+1)]
                #sample2_coverage = coverage.loc[coverage["Id"] == id2]
                #sample2_coverage = sample2_coverage.loc[sample2_coverage["chr"] == "Sabin2ref"]
                #sample2_coverage = sample2_coverage.loc[sample2_coverage["chr.pos"] == str(i+1)]
                #cov1 = int(sample1_coverage["coverage"])
                #cov2 = int(sample2_coverage["coverage"])
                #print("Warning! Bases at position " + str(i+1) + " in " + id1 + " and " + id2 + " are different! Using the replicate with higher coverage at the site: " + str(cov1) + "x vs. " + str(cov2) + "x")
                #if(cov1 >= cov2):
                #    new_base = base1
                #elif(cov2 > cov1):
                #    new_base = base2
                #else:
                #    print("Something went terribly wrong! Abort!")

        new_sequence[i] = new_base

    new_sequence_string = ''.join(new_sequence)

    return new_sequence_string

def PrintNewFasta(resolved, bio_sample_id, outdir):

    new_filename = outdir + bio_sample_id + ".removed.parsed.filtered.resolved.fasta"
    #new_filename = outdir + bio_sample_id + ".filtered.resolved.fasta"
    record = SeqRecord(Seq(resolved), id = bio_sample_id)
    SeqIO.write(record, new_filename, "fasta")

    return

# ========================= Main =============================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', action="store", dest="indir")
    parser.add_argument('--outdir', action="store", dest="outdir")
    parser.add_argument('--meta', action="store", dest="metafile")
    parser.add_argument('--coverage', action="store", dest="coverage")
    args = parser.parse_args()

    metadata = pd.read_csv(args.metafile, index_col = None, header = 0, dtype = object)
    metadata.columns = metadata.columns.astype(str)

    coverage = LoadCoverage(args.coverage)

    for index, row in metadata.iterrows():
        if(int(row["reps"]) == 2):
            id1 = str(row["sequencingID_1"])
            id2 = str(row["sequencingID_2"])
            bio_sample_id = str(row["anonSpecId"])
            file1 = args.indir + id1 + ".removed.parsed.filtered.fasta" # deepSNV version
            file2 = args.indir + id2 + ".removed.parsed.filtered.fasta"
            #file1 = args.indir + id1 + ".filtered.fasta" # manual version
            #file2 = args.indir + id2 + ".filtered.fasta"
            file1path = Path(file1)
            file2path = Path(file2)
            #print("Working on duplicates for biological sample: " + bio_sample_id)
            if(file1path.is_file() is False or file2path.is_file() is False):
                print("We don't have both replicates.")
                continue
            seq1 = SeqIO.parse(file1, "fasta")
            seq2 = SeqIO.parse(file2, "fasta")
            resolved = ResolveFasta(seq1, seq2, coverage, id1, id2, bio_sample_id)
            if(resolved is None):
                print("No resolved seq for " + bio_sample_id)
                continue
            #PrintNewFasta(resolved, bio_sample_id, args.outdir)
        elif(int(row["reps"]) == 1):
            id = str(row["sequencingID_1"])
            bio_sample_id = str(row["anonSpecId"])
            infile = args.indir + id + ".removed.parsed.filtered.fasta"
            outfile = args.outdir + bio_sample_id + ".removed.parsed.filtered.resolved.fasta"
            #infile = args.indir + id + ".filtered.fasta"
            #outfile = args.outdir + bio_sample_id + ".filtered.resolved.fasta"
            command = "cp " + infile + " " + outfile
            #os.system(command)
        else:
            print("Something is wrong; reps should be 1 or 2.")
            break


if __name__ == "__main__":
    main()
