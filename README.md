# DMFlux

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-20xx.xxxxx%20-green.svg)](https://arxiv.org/abs/20xx.xxxxx)

DMFlux is a package to organize calculations of neutrinos from dark matter annihilation/decay.

![DMFlux](docs/source/logo.jpg)

## Overview
Indirect detection which detects Standard Model (SM) particles produced by dark matter annihilation/decay is an important piece of current approaches searching for the dark matter. Stable particles from astrophysical sources are messengers of these indirect signals. Among messengers used in indirect searches for dark matter, neutrinos are special as they are neutral, light, and seldom interact. These unique properties give them advantages in astrophysical studies: they are advantageous over cosmic rays as they can point back to their sources and unlike gamma rays can exit environments of large matter and radiation densities. It is important to have a tool to organize the neutrino flux generation.

## Features
* **DM Capture**: Capture and annihilation rate calculation of dark matter.  
	Things related to dark matter capture are in the [`capture`](capture/) folder.
* **Flux at Production**: We provide neutrino fluxes at production using PYTHIA8.2 with electroweak correction from [`Bauer et al.`](https://arxiv.org/abs/20xx.xxxxx).
    Things related to neutrino fluxes at productions are in the [`production`](production/) folder. Fluxes with and without electroweak correction are in [`production/data`](production/data/).
* **Propagation**: Propagate initial fluxs to Earth or the detector.     
	Things related to neutrino propagation are in the [`propagate`](propagate/) folder.
* **Secluded DM**: Besides the standard case, DMFlux also includes the possibility of a secluded dark matter sector which introduces a long-lived mediator. 
 

## Dependencies
* DMFlux needs a modified version of nuSQuIDs which can be installed from 
  [`nuSQuIDS branch`](https://github.com/qrliu/nuSQuIDS) and please make sure python interface is activated.
* Python
* NumPy
* Scipy
* SymPy



## Installation
DMFlux can be installed using the following command
```
python setup.py install
```

If you want to generate flux at production yourself, please go to [`production`](production/) and see [`README`](production/README.md) file there. 
 
## Examples
You can find examples of how to run DMFlux in the 
[`examples`](examples) folder.
