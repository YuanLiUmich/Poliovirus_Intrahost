mkdir Benchmark
mkdir Benchmark/data
mkdir Benchmark/data/fastq
mkdir Benchmark/data/reference
# copy over reference files here
cp /nfs/turbo/med-alauring2/raw_data/2019/Run_2994/lauring/*/* data/fastq/
python change_miseq_names_Matlab.py -s data/fastq -f data/fastq_renamed -run
cd data/fastq_renamed
gunzip -v *gz
mv 133396.1.1.fastq Sabin1_PlasmidControl.1.1.fastq
mv 133396.2.1.fastq Sabin1_PlasmidControl.2.1.fastq
cd ../../
# Run alignment stage to concat combined reference: pipeline_align.sbat
# Subset aligned reads with cap_bam.txt
# Then run deepSNV stage on capped reads with pipeline_capped.sbat
