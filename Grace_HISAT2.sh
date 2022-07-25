#!/bin/bash
#SBATCH --export=NONE               				# do not export current env to the job
#SBATCH --job-name=hisat2           				# job name
#SBATCH --time=5-00:00:00           				# max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         				# tasks (commands) per compute node
#SBATCH --cpus-per-task=28          				# CPUs (threads) per command
#SBATCH --mem=56G                   				# total memory per node
#SBATCH --output=stdout.%j          				# save stdout to file
#SBATCH --error=stderr.%j           				# save stderr to file
#SBATCH --mail-type=ALL              				#Send email on all job events
#SBATCH --mail-user=venura.herath@tamu.edu 			#Email Notification

#Grace

module load GCC/9.3.0 OpenMPI/4.0.3 HISAT2/2.2.1


########## INPUTS ##########
Mock="TM_MOCK_2hr_1 TM_MOCK_2hr_2 TM_MOCK_2hr_3 TM_MOCK_5hr_1 TM_MOCK_5hr_2 TM_MOCK_5hr_3"
TM="TM_2hr_1 TM_2hr_2 TM_2hr_3 TM_5hr_1 TM_5hr_2 TM_5hr_3" 

# you can use an already prefixed genome found at: /scratch/data/bio/genome_indexes/
genome_index_prefix='/scratch/data/bio/genome_indexes/other_genomes/potato/castle_russet/hisat2/CR_v2.0_pseudomolecules'

######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK


############Commands#################

for Mock in $Mock; do
    hisat2 -p $threads --rg-id=${Mock} --rg SM:Mock --dta --rna-strandness RF -x $genome_index_prefix -1 ${Mock}_1.fq.gz -2 ${Mock}_2.fq.gz -S ./RNA_ALIGN/${Mock}.sam
done

for TM in $TM; do
    hisat2 -p $threads --rg-id=${TM} --rg SM:TM --dta --rna-strandness RF -x $genome_index_prefix -1 ${TM}_1.fq.gz -2 ${TM}_2.fq.gz -S ./RNA_ALIGN/${TM}.sam
done


