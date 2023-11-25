#!/bin/bash -l

#SBATCH --job-name=job_name
#SBATCH --partition=job_partition
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=1G
#SBATCH --output=%x_%j.output
#SBATCH --error=%x_%j.error
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email

source /home/${USER}/.bashrc
conda activate pbilby #activate your conda environment

fname=xibar_GW170817_binary_love.py #parameter estimation script
strain_GW170817=$HOME/scratch/GW170817-Data/no-glitch #path to the data

#switch to scratch directory and copy over files from your project directory
savedir=$HOME/scratch/$SLURM_JOB_NAME/$SLURM_JOB_ID
olddir=$PWD
mkdir -p $savedir
cd  $savedir
cp $olddir/injection.py .
cp $olddir/$fname .
cp $olddir/nrtidal_d.py .
cp $olddir/updated_binary_love_marginalized.py .
cp $olddir/sampler_params.py .

python $fname -n $SLURM_CPUS_PER_TASK -sd $strain_GW170817

