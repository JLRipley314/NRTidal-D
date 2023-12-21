#!/bin/bash -l

#
# Usage:
#
# sbatch launch.slurm [value of xitilde] 
#

#SBATCH --job-name=GW17

# #SBATCH --partition=GravityTheory
# #SBATCH --partition=caps

#SBATCH --partition=IllinoisComputes
#SBATCH --account=nyunes-ic

#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=1G
#SBATCH --output=%x_%j.output
#SBATCH --error=%x_%j.error

# #SBATCH --mail-user=ripley@illinois.edu
# #SBATCH --mail-type=FAIL,END

export OMP_NUM_THREADS=1
export MPI_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Assuming you have a conda environment called "bilby"
# See the environment.yml file
source /home/${USER}/.bashrc
source activate bilby 

python main.py -n $SLURM_CPUS_PER_TASK -x ${1}
