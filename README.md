# Analysis of within-host diversity of OPV2

This repository holds code and data for the analysis and relies heavily on our other repository, [variant_pipeline](https://github.com/lauringlab/variant_pipeline), for primary variant calling. The rest of the analysis scripts are provided here.

# Overview
--------

    project
    |- README          # The top level description of content. You are here.
    |
    |- data  
    |  |- metadata/  # Sample metadata for the specimens sequenced in this study.
    |  |- reference/  # Reference fasta files used in sequence read alignment, accessions used for comparison to cVDPVs, and other reference info.
    |  |- raw/         # Raw data, organized by the eight sequencing runs (96-well plates). Includes coverage data, raw consensus files, and raw variant calls relative to sample consensus. Variant calls with coding changes denoted relative to the reference are in plate*ctl/Final_variants.
    |  |- processed/     # Cleaned and intermediate data. Please note that the variant files are large (> 100 MB) and are not uploaded to this repository.
    |  |- benchmarking/     # Data and scripts for benchmarking analysis of variant calling from defined mixtures.
    |- scripts/           # All analysis code for trial specimens.
    |  |- primary_analysis/    # The HPC job scripts and option files used to process the sequencing data from fastq format to variant calls. Each plate has a file describing the steps from raw data through running the pipeline.
    |  |- secondary_analysis/  # R and python scripts used for variant analysis after primary variant calling.
    
  --------

# Notes

Raw sequencing data is available through the NCBI SRA at BioProject PRJNA637613.

The [variant pipeline](https://github.com/lauringlab/variant_pipeline) was used for sequence alignment and variant calling with deepSNV. 
The adapter trimming and read mapping options were modified in a local version of [this configuration file](https://github.com/lauringlab/variant_pipeline/blob/master/scripts/variantPipeline.bpipe.stages.groovy).

# Contact

If you have questions, please contact the [Lauring Lab](https://lauringlab.wordpress.com/contacts/).