#!/bin/bash
#SBATCH -J gif #job name
#SBATCH --time=1:00:00 #requested time (DD-HH:MM:SS)
#SBATCH -p batch #running on "mpi" partition/queue
#SBATCH -N 1 #1 nodes
#SBATCH -n 1 #1 tasks total (default 1 CPU core per task) = # of cores
#SBATCH --mem=4g #requesting 2GB of RAM total
#SBATCH --output=logs/gif_%a.%j.%N.out
#SBATCH --error=logs/gif_%a.%j.%N.err
#SBATCH --mail-type=ALL #email options
#SBATCH --mail-user=wwhite06@tufts.edu
 
# get line number ${SLURM_ARRAY_TASK_ID} from tasks file
CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" gif_cmds.sh)
# tell bash to run $CMD
echo "${CMD}" | bash
