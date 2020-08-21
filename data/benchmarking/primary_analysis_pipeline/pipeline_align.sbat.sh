#!/bin/bash

######### Slurm Preamble #########################################

#“#SBATCH” directives that convey submission options:
##### The name of the job
#SBATCH --job-name=cap

##### When to send e-mail: pick from NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=BEGIN,END,FAIL

##### Resources for your job: number of physical nodes
#SBATCH --nodes=1

##### Resources for your job: number of task per a node (number of CPU-cores per a node)
#SBATCH --ntasks-per-node=2

##### Resources for your job: memory per a CPU-core
#SBATCH --mem-per-cpu=32gb

##### Maximum amount of time the job will be allowed to run. Recommended formats: MM:SS, HH:MM:SS, DD-HH:MM
#SBATCH --time=8:00:00

##### The resource account; who pays
#SBATCH --account=alauring

##### Output Name
#SBATCH --output=pipeline_capped.log

########## End of preamble! #########################################
# No need to “cd”. Slurm starts the job in the submission directory.
#####################################################################
# The application(s) to execute along with its input arguments and options:

/bin/hostname

rm -r /home/avalesan/.bpipedb/jobs/
python ~/variant_pipeline/bin/variantPipeline.py ~/variant_pipeline/scripts/aligning_pipeline.groovy "./data/fastq_renamed/*fastq" ./data/aligned_output ./options_shiny_onesided.yaml

echo done!
echo Finish time: `date`
