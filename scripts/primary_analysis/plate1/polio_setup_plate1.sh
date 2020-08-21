cd /scratch/alauring_fluxm/avalesan/
mkdir Matlab_OPV2
cd Matlab_OPV2
mkdir plate1_Run_2894
cd plate1_Run_2894
mkdir data
mkdir data/fastq
mkdir data/reference
# copy over or upload reference files
cp /nfs/turbo/med-alauring2/raw_data/2019/Run_2894/lauring/*/* data/fastq/
python ../change_miseq_names_Matlab.py -s data/fastq -f data/fastq_renamed -run
cd data/fastq_renamed
gunzip -v *gz
mv 131650.1.1.fastq Sabin2_PlasmidControl.1.1.fastq
mv 131650.2.1.fastq Sabin2_PlasmidControl.2.1.fastq
cd ../../
qsub submit_pipeline.pbs
