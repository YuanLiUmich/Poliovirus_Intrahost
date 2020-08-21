echo 2A
muscle -in MatlabBioSampleConsensus_mOPV2_2A.fasta -physout 2A.phylip
cat 2A.phylip | sed 's/P........./&\n/g' > 2A.space.phylip

echo 2B
muscle -in MatlabBioSampleConsensus_mOPV2_2B.fasta -physout 2B.phylip
cat 2B.phylip | sed 's/P........./&\n/g' > 2B.space.phylip

echo 2C
muscle -in MatlabBioSampleConsensus_mOPV2_2C.fasta -physout 2C.phylip
cat 2C.phylip | sed 's/P........./&\n/g' > 2C.space.phylip

echo 3A
muscle -in MatlabBioSampleConsensus_mOPV2_3A.fasta -physout 3A.phylip
cat 3A.phylip | sed 's/P........./&\n/g' > 3A.space.phylip

echo 3B
muscle -in MatlabBioSampleConsensus_mOPV2_3B.fasta -physout 3B.phylip
cat 3B.phylip | sed 's/P........./&\n/g' > 3B.space.phylip

echo 3C
muscle -in MatlabBioSampleConsensus_mOPV2_3C.fasta -physout 3C.phylip
cat 3C.phylip | sed 's/P........./&\n/g' > 3C.space.phylip

echo 3D
muscle -in MatlabBioSampleConsensus_mOPV2_3D.fasta -physout 3D.phylip
cat 3D.phylip | sed 's/P........./&\n/g' > 3D.space.phylip

echo VP1
muscle -in MatlabBioSampleConsensus_mOPV2_VP1.fasta -physout VP1.phylip
cat VP1.phylip | sed 's/P........./&\n/g' > VP1.space.phylip

echo VP2
muscle -in MatlabBioSampleConsensus_mOPV2_VP2.fasta -physout VP2.phylip
cat VP2.phylip | sed 's/P........./&\n/g' > VP2.space.phylip

echo VP3
muscle -in MatlabBioSampleConsensus_mOPV2_VP3.fasta -physout VP3.phylip
cat VP3.phylip | sed 's/P........./&\n/g' > VP3.space.phylip

echo VP4
muscle -in MatlabBioSampleConsensus_mOPV2_VP4.fasta -physout VP4.phylip
cat VP4.phylip | sed 's/P........./&\n/g' > VP4.space.phylip

echo done!
