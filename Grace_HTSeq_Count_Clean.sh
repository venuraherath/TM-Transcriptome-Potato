#!/bin/bash
#SBATCH --export=NONE               								# do not export current env to the job
#SBATCH --job-name=HTSEQ_Count           							# job name
#SBATCH --time=4-00:00:00           								# max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         								# tasks (commands) per compute node
#SBATCH --cpus-per-task=28          								# CPUs (threads) per command
#SBATCH --mem=256G                   								# total memory per node
#SBATCH --output=stdout.%j          								# save stdout to file
#SBATCH --error=stderr.%j           								# save stderr to file
#SBATCH --mail-type=ALL              								#Send email on all job events
#SBATCH --mail-user=venura.herath@tamu.edu 							#Email Notification

#Grace
module load GCC/11.2.0  OpenMPI/4.1.1 HTSeq/2.0.1

########## INPUTS ##########
Mock="TM_MOCK_2hr_1 TM_MOCK_2hr_2 TM_MOCK_2hr_3 TM_MOCK_5hr_1 TM_MOCK_5hr_2 TM_MOCK_5hr_3"
TM="TM_2hr_1 TM_2hr_2 TM_2hr_3 TM_5hr_1 TM_5hr_2 TM_5hr_3"


######## PARAMETERS ########
#threads=$SLURM_CPUS_PER_TASK

################CODE################

mkdir htseq_counts_clean

for Mock in $Mock; do
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id ${Mock}.bam cr.working_models.pm.locus_assign_clean.gtf > htseq_counts_clean/${Mock}.tsv
done

for TM in $TM; do
	htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id ${TM}.bam cr.working_models.pm.locus_assign_clean.gtf > htseq_counts_clean/${TM}.tsv
done

