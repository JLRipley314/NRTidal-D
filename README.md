# NRTidal-D

Analysis of the dissipative tidal deformability from gravitational wave (GW) strain data. 
In short: scripts for the Bayesian parameter estimation of the dissipative tidal deformabilities of the neutron stars in GW170817.
For more details of our analysis, see the accompanying paper (see the `Citation` at the bottom). 

1. `Waveform model`: Scripts for our waveform model (including the marginalized binary Love relations) 
2. `Injection-Recovery`: Scripts to run Bayesian analysis of injection/recovery analysis for different GW detector networks.  
3. `GW170817-Recovery`: Scripts to perform Bayesian analysis of GW170817 GW strain data. 
3. `ASD-Files`: Amplitude strain data files (for different gravitational wave detector networks). 

# Installation

This is not a proper python package; it is a collection of scripts we used to perform Bayesian analysis on GW data. 
To run these scripts, we recommend making a new `conda` environment. 

The `environment.yml` file contains more details about the dependencies of our scripts.
If you have any problems getting the environment set up, please email us.
The main dependency is the [bilby](https://lscsoft.docs.ligo.org/bilby/) gravitational wave analysis library. 
If you have that installed our scripts should run fine, but beware of the different versions of bilby; not everything may work properly if you use an incompatible version.

# GW170817 Data 

We used the glitch removed strain data, which can be accessed [here](https://gwosc.org/events/GW170817/).

# Authors

If you have any questions, please feel free to contact any one of us. 

[Justin Ripley](https://github.com/JLRipley314)

[Abhishek Hegade](https://github.com/AbhiHegade) 

[Rohit Chandramouli](https://github.com/rsc496)

# Citation

We will add the appropriate citation once the paper appears on the arxiv.
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
