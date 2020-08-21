cd /scratch/alauring_fluxm/avalesan/Matlab_OPV2
mkdir plate2_Run_2905
cd plate2_Run_2905
mkdir data
mkdir data/fastq
mkdir data/reference
# copy over reference files
cp /nfs/turbo/med-alauring2/raw_data/2019/Run_2905/lauring/*/* data/fastq/
python ../change_miseq_names_Matlab.py -s data/fastq -f data/fastq_renamed -run
cd data/fastq_renamed
gunzip -v *gz
mv 131747.1.1.fastq Sabin2_PlasmidControl.1.1.fastq
mv 131747.2.1.fastq Sabin2_PlasmidControl.2.1.fastq
cd ../../
qsub submit_pipeline.pbs
