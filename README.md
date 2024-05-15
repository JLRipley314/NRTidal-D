# NRTidal-D

Analysis of the dissipative tidal deformability from gravitational wave (GW) strain data. 
In short: scripts for the Bayesian parameter estimation of the dissipative tidal deformabilities of the neutron stars in GW170817.
For more details of our analysis, see the accompanying paper (see the `Citation` at the bottom). 

1. `Waveform model`: Scripts for our waveform model (including the marginalized binary Love relations) 
2. `Injection-Recovery`: Scripts to run Bayesian analysis of injection/recovery analysis for different GW detector networks.  
3. `GW170817-Recovery`: Scripts to perform Bayesian analysis of GW170817 GW strain data. 
3. `ASD-Files`: Amplitude strain data files (for different gravitational wave detector networks). 

# Conventions

In our ``theory'' and ``analysis'' papers (see Citation), we define the dissipative tidal deformability ``xi_tilde'' with an extra factor of 8, so that in the limit where the stars have the same radii, masses, tidal deformability, and dissipative tidal deformability (R1 = R2, ..., Xi1 = Xi2), then Xi reduces to Xi1. 
In the code here, we did not include the extra factor of 8 for historical reasons.
To match the output in our papers then, you should multiply the code output for xi_tilde by 8. 

# Installation

This is not a proper python package; it is a collection of scripts we used to perform Bayesian analysis on GW data. 
To run these scripts, we recommend making a new `conda` environment. 

The `environment.yml` file contains more details about the dependencies of our scripts.
If you have any problems getting the environment set up, please email us.
The main dependency is the [bilby](https://lscsoft.docs.ligo.org/bilby/) gravitational wave analysis library. 
If you have that installed our scripts should run fine, but beware of the different versions of bilby; not everything may work properly if you use an incompatible version.

# GW170817 Data 

We downloaded the glitch-free GW170817 strain data from [gwosc.org](https://doi.org/10.7935/K5B8566F).

# Data, plots 

We have uploaded the strain data we used, the raw output of our samplering of that data, and our plotting scripts to [zenodo](https://doi.org/10.5281/zenodo.10967278).

# Authors

If you have any questions, please feel free to contact any one of us. 

[Justin Ripley](https://github.com/JLRipley314)

[Abhishek Hegade](https://github.com/AbhiHegade) 

[Rohit Chandramouli](https://github.com/rsc496)

# Citation

Our data analysis paper 
```
@article{Ripley:2023lsq,
    author = "Ripley, Justin L. and Hegade K. R., Abhishek and Chandramouli, Rohit S. and Yunes, and Nicolas",
    title = "{First constraint on the dissipative tidal deformability of neutron stars}",
    eprint = "2312.11659",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "12",
    year = "2023"
}
```

For reference, the "theory" paper (derivation of the dissiaptive tidal deformability phasing formula) is
```
@article{Ripley:2023qxo,
    author = "Ripley, Justin L. and Hegade K. R., Abhishek and Yunes, Nicolas",
    title = "{Probing internal dissipative processes of neutron stars with gravitational waves during the inspiral of neutron star binaries}",
    eprint = "2306.15633",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    doi = "10.1103/PhysRevD.108.103037",
    journal = "Phys. Rev. D",
    volume = "108",
    number = "10",
    pages = "103037",
    year = "2023"
}
``` 
