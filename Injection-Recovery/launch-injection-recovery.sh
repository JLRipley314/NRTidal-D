#!/usr/bin/bash

#
# Launches multiple injection/recovery runs. 
# Saves to outstem.
#

nrtd=$HOME/NRTidal-D

#
# Manually set this
#
outstem=/home/ripley/scratch/GW170817-InjectionRecovery

detectors=("O4" "O5" "CE")
xitildes=(20 200 400)

for det in "${detectors[@]}"; do
  for xitilde in "${xitildes[@]}"; do
    echo "Detector ${det}, Xitilde ${xitilde}"
    main="main_${det}.py"
    d="$outstem/${det}-recover-Xi${xitilde}"
    mkdir $d
    cd $d
    cp $nrtd/Waveform-Model/nrtidal_d.py nrtidal_d.py 
    cp $nrtd/Waveform-Model/updated_binary_love_marginalized.py updated_binary_love_marginalized.py 
    cp $nrtd/Injection-Recovery/$main main.py 
    cp $nrtd/ASD-Files/*.txt . 
    cp $nrtd/Injection-Recovery/launch.slurm launch.slurm
    sbatch launch.slurm ${xitilde} 
  done
done
