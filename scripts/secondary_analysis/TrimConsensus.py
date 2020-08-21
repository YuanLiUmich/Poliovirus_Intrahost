
### Author: Andrew Valesano
### Project: Matlab Poliovirus Sequencing
### Purpose: Use coverage information to trim each consensus sequence from the pipeline.
### Only returns a file if all of the positions between the ends are also above the cutoff.

### Usage: python ./scripts/TrimConsensus.py --indir data/raw/ --outdir data/processed/consensus_trimmed/ --cutoff 20 --genome Sabin2ref

# ======================= Import modules ======================

import argparse
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ========================= Functions =============================

def LoadCoverage(indir):

    coverage_files = "%s*/all.coverage.csv" % indir
    cov = []
    for file in glob.glob(coverage_files):
        print(file)
        c = pd.read_csv(file, index_col = None, header = 0)
        cov.append(c)
    coverage = pd.concat(cov, axis=0, ignore_index=True)
    coverage.columns = coverage.columns.astype(str)

    return coverage

def TrimFasta(filename, coverage, cutoff, genome):

    # Get coverage for the sample
    filename_only = filename.split("/")[4]
    id = filename_only.split(".")[0]
    sample_coverage = coverage.loc[coverage["Id"] == id]
    sample_coverage = sample_coverage.loc[sample_coverage["chr"] == genome]

    # Find start and end positions where coverage is at least at the cutoff
    coverage_above_cutoff = sample_coverage.loc[sample_coverage["coverage"] >= cutoff]
    if coverage_above_cutoff.empty:
        return None
    cov_values = coverage_above_cutoff["coverage"]
    if not all(num > cutoff for num in cov_values):
        return None
    start = min(coverage_above_cutoff["concat.pos"])
    end = max(coverage_above_cutoff["concat.pos"])

    # Trim ends based on coverage
    fastaFile = SeqIO.parse(filename, "fasta")
    for entry in fastaFile:
        if entry.id == genome:
            trimmed = entry.seq[start:end]
            trimmedSeq = SeqRecord(trimmed, entry.id)

    return trimmedSeq

def PrintNewFasta(trimmed, filename, outdir):

    filename_only = filename.split("/")[4]
    id = filename_only.split(".")[0]
    new_filename = outdir + id + ".removed.parsed.trimmed.fasta"
    SeqIO.write(trimmed, new_filename, "fasta")

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
    for file in glob.glob(infiles):
        print("Working on file: " + file)
        trimmed = TrimFasta(file, coverage, args.cutoff, args.genome)
        if trimmed is None:
            continue
        PrintNewFasta(trimmed, file, args.outdir)


if __name__ == "__main__":
    main()
