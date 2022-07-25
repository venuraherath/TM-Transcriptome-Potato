#!/bin/bash
#SBATCH --export=NONE               				# do not export current env to the job
#SBATCH --job-name=hisat2_samtools           				# job name
#SBATCH --time=5-00:00:00           				# max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         				# tasks (commands) per compute node
#SBATCH --cpus-per-task=28          				# CPUs (threads) per command
#SBATCH --mem=56G                   				# total memory per node
#SBATCH --output=stdout.%j          				# save stdout to file
#SBATCH --error=stderr.%j           				# save stderr to file
#SBATCH --mail-type=ALL              				#Send email on all job events
#SBATCH --mail-user=venura.herath@tamu.edu 			#Email Notification

module load GCC/8.3.0 GCC/9.3.0 iccifort/2019.5.281 SAMtools/1.10

########## INPUTS ##########
Mock="TM_MOCK_2hr_1 TM_MOCK_2hr_2 TM_MOCK_2hr_3 TM_MOCK_5hr_1 TM_MOCK_5hr_2 TM_MOCK_5hr_3"
TM="TM_2hr_1 TM_2hr_2 TM_2hr_3 TM_5hr_1 TM_5hr_2 TM_5hr_3" 

######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK

########################SAMTOOLS######################

############Commands#################
cd RNA_ALIGN

for Mock in $Mock; do
    samtools sort -@ $threads -o ${Mock}.bam ${Mock}.sam
done

for TM in $TM; do
    samtools sort -@ $threads -o ${TM}.bam ${TM}.sam
done




