cd /scratch/alauring_fluxm/avalesan/Matlab_OPV2
mkdir plate5_Run_2946
cd plate5_Run_2946
mkdir data
mkdir data/fastq
mkdir data/reference
# copy over reference files here
cp /nfs/turbo/med-alauring2/raw_data/2019/Run_2946/lauring/*/* data/fastq/
python ../change_miseq_names_Matlab.py -s data/fastq -f data/fastq_renamed -run
cd data/fastq_renamed
gunzip -v *gz
mv 132038.1.1.fastq Sabin2_PlasmidControl.1.1.fastq
mv 132038.2.1.fastq Sabin2_PlasmidControl.2.1.fastq
cd ../../
qsub submit_pipeline.pbs
