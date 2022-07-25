#!/bin/bash                                                        
#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=fastqc           # job name       
#SBATCH --time=10:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=2           # CPUs (threads) per command  
#SBATCH --mem=4G                    # total memory per node
#SBATCH --output=stdout.%j          # save stdout to file
#SBATCH --error=stderr.%j           # save stderr to file
#SBATCH --mail-type=ALL             #Send email on all job events
#SBATCH --mail-user=venura.herath@tamu.edu #Email Notification

######Grace#########

module load FastQC/0.11.9-Java-11

########## INPUTS ##########
SAMPLES="TM_2hr_1 TM_2hr_2 TM_2hr_3 TM_5hr_1 TM_5hr_2 TM_5hr_3 TM_MOCK_2hr_1 TM_MOCK_2hr_2 TM_MOCK_2hr_3 TM_MOCK_5hr_1 TM_MOCK_5hr_2 TM_MOCK_5hr_3"

######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK

################################### COMMANDS ###################################
mkdir fastQC
for SAMPLE in $SAMPLES; do
fastqc -t $threads -o ./fastQC ${SAMPLE}_1.fq.gz ${SAMPLE}_2.fq.gz
done
################################################################################