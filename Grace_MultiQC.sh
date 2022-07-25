#!/bin/bash                                                        
#SBATCH --export=NONE               		# do not export current env to the job
#SBATCH --job-name=multiqc           		# job name       
#SBATCH --time=05:00:00             		# max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         		# tasks (commands) per compute node
#SBATCH --cpus-per-task=2           		# CPUs (threads) per command  
#SBATCH --mem=4G                    		# total memory per node
#SBATCH --output=stdout.%j          		# save stdout to file
#SBATCH --error=stderr.%j           		# save stderr to file
#SBATCH --mail-user=venura.herath@tamu.edu 	#Email Notification


#Grace

module load GCC/9.3.0 OpenMPI/4.0.3 iccifort/2020.1.217 impi/2019.7.217 MultiQC/1.9-Python-3.8.2

######## PARAMETERS ########
threads=$SLURM_CPUS_PER_TASK

#######################Commands##############
python3 -m multiqc .


