#!/usr/bin/bash

nrtd=$HOME/NRTidal-D

#
# Manually set this
#
outstem=/home/ripley/scratch/GW170817-InjectionRecovery

DETECTORS=("O4" "O5" "CE")
XIBARS=(20 200 400)

for det in "${DETECTORS[@]}"; do
  for xibar in "${XIBARS[@]}"; do
    echo "Detector ${det}, Xibar ${xibar}"
    main="main_${det}.py"
    d="$outstem/${det}-recover-Xi${xibar}"
    mkdir $d
    cd $d
    cp $nrtd/Waveform-Model/nrtidal_d.py nrtidal_d.py 
    cp $nrtd/Waveform-Model/updated_binary_love_marginalized.py updated_binary_love_marginalized.py 
    cp $nrtd/Injection-Recovery/$main main.py 
    cp $nrtd/ASD-Files/*.txt . 
    cp $nrtd/Injection-Recovery/launch.slurm launch.slurm
    sbatch launch.slurm ${xibar} 
  done
done
