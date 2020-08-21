
### Author: Andrew Valesano
### Project: Matlab Poliovirus Sequencing
### Purpose: Use coverage information to filter each consensus sequence from the pipeline (parsed_fa).
### If a position in the fasta has coverage below the cutoff, an N will be placed there.
### Working directory: Poliovirus_Intrahost

### Usage: python ./scripts/secondary_analysis/FilterConsensus.py --indir data/raw_final20/ --outdir data/processed/consensus_filtered_final/ --cutoff 10 --genome Sabin2ref

# ======================= Import modules ======================

import argparse
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

def FilterFasta(filename, coverage, cutoff, genome):

    # Get coverage for the sample
    #filename_only = filename.split("/")[2] # 4 for deepSNV version
    filename_only = filename.split("/")[4] # 4 for deepSNV version
    id = filename_only.split(".")[0]
    sample_coverage = coverage.loc[coverage["Id"] == id]
    sample_coverage = sample_coverage.loc[sample_coverage["chr"] == genome]

    if(len(sample_coverage) == 0):
        print("Can't find coverage data for this one.")
        return None

    # Loop through original fasta, building new sequence string.
    # Coverage is from 1 to 7439, while entry.seq starts at 0.
    fastaFile = SeqIO.parse(filename, "fasta")
    for entry in fastaFile:
        if entry.id == genome:
            genomeLength = len(entry.seq)
            new_sequence = list("x" * genomeLength)
            for i in range(1, genomeLength+1, 1):
                cov_at_position = sample_coverage.loc[sample_coverage["chr.pos"] == str(i)]
                new_base = "x"
                if(len(cov_at_position) == 0):
                    new_base = "N"
                elif(int(cov_at_position["coverage"]) < cutoff):
                    new_base = "N"
                else:
                    new_base = entry.seq[i-1]
                new_sequence[i-1] = new_base
            new_sequence_string = ''.join(new_sequence)

    return new_sequence_string

def PrintNewFasta(filtered, filename, outdir, genome, cutoff):

    #filename_only = filename.split("/")[2] # 2 for samtools version
    filename_only = filename.split("/")[4] # 4 for deepSNV version
    id = filename_only.split(".")[0]
    new_filename = outdir + id + ".removed.parsed.filtered.fasta"
    #new_filename = outdir + id + ".filtered.fasta"
    record = SeqRecord(Seq(filtered), id = id)
    SeqIO.write(record, new_filename, "fasta")

    return

# ========================= Main =============================

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', action="store", dest="indir")
    parser.add_argument('--outdir', action="store", dest="outdir")
    parser.add_argument('--cutoff', action="store", dest="cutoff", type=int)
    parser.add_argument('--genome', action="store", dest="genome", default = "Sabin2ref")
    args = parser.parse_args()

    coverage = LoadCoverage(args.indir)

    infiles = "%s*/parsed_fa/*.fasta" % args.indir
    #infiles = "data/consensus/*.fa"
    count = 1
    for file in glob.glob(infiles):
        filepath = Path(file)
        if(filepath.is_file() is False):
            count = count + 1
            print("File " + file + " isn't there.")
            continue
        #filename_only = file.split("/")[2] # 4 for deepSNV version
        filename_only = file.split("/")[4] # 4 for deepSNV version
        id = filename_only.split(".")[0]
        if("Sabin2_PlasmidControl" in id):
            print("Skipping the plasmid control.")
            continue

        #new_filename = args.outdir + id + ".filtered.fasta"
        new_filename = args.outdir + id + ".removed.parsed.filtered.fasta"
        newfilepath = Path(new_filename)
        if(newfilepath.is_file() is True):
            count = count + 1
            print("File " + file + " is already done.")
            continue

        print("Working on file: " + file + ", file number " + str(count))
        filtered = FilterFasta(file, coverage, args.cutoff, args.genome)
        if(filtered is None):
            print("File " + file + " is empty ")
            count = count + 1
            continue
        PrintNewFasta(filtered, file, args.outdir, args.genome, args.cutoff)
        count = count + 1
    print("Final count: " + str(count))

if __name__ == "__main__":
    main()
