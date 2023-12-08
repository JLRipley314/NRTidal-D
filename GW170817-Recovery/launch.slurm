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

nrtd=$HOME/NRTidal-D

source /home/${USER}/.bashrc
conda activate bilby #activate your conda environment

#
# Manually set this: path to the strain data
#
strain_GW170817=$HOME/scratch/GW170817-Data/no-glitch 

#
# Manually set this: path to directory to save to
#
savedir=$HOME/scratch/$SLURM_JOB_NAME/$SLURM_JOB_ID

# Parameter estimation script
fname=xibar_GW170817_binary_love_IMRPhenomPv2.py 
#fname=xibar_GW170817_binary_love_IMRPhenomD.py 
#fname=xibar_GW170817_no_binary_love.py 

mkdir -p $savedir
cd  $savedir
cp $nrtd/GW170817-Recovery/$fname .
cp $nrtd/Waveform-Model/nrtidal_d.py .
cp $nrtd/Waveform-Model/updated_binary_love_marginalized.py .

python $fname -n $SLURM_CPUS_PER_TASK -sd $strain_GW170817
