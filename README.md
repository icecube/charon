# χarον

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-20xx.xxxxx%20-green.svg)](https://arxiv.org/abs/20xx.xxxxx)

χarον is a package to organize calculations of neutrinos from dark matter annihilation/decay.

<!---[χαρον](docs/source/logo.jpg)-->

## Overview
Indirect detection which detects Standard Model (SM) particles produced by dark matter annihilation/decay is an important piece of current approaches searching for the dark matter. Stable particles from astrophysical sources are messengers of these indirect signals. Among messengers used in indirect searches for dark matter, neutrinos are special as they are neutral, light, and seldom interact. These unique properties give them advantages in astrophysical studies: they are advantageous over cosmic rays as they can point back to their sources and unlike gamma rays can exit environments of large matter and radiation densities. It is important to have a tool to organize the neutrino flux generation.

## Features
* **Neutrino Flux at Production**: We provide neutrino fluxes at production using PYTHIA8.2 with electroweak corrections from [`BRW`](https://arxiv.org/abs/20xx.xxxxx) for DM mass larger than 500 GeV and the option to generate fluxes without EW corrections onself.
    Things related to neutrino fluxes at productions are in the [`production`](production/) folder.  
* **Neutrino Propagation**: Propagate initial fluxs to different locations. Tables of initial neutrino fluxes with and without electroweak correction are in [`charon/data`](charon/data/). When using χarον to propagate, the default is that when DM mass is smaller than 500 GeV the initial flux is without EW correction and when DM mass is larger than 500 GeV the initial flux has EW correction included.      
* **DM Capture**: Some models of DM capture and annihilation rate calculation.  
* **J Factor**: Compute J factor.
* **Secluded DM**: Besides the standard case, χarον also includes the possibility of a secluded dark matter sector which introduces a long-lived mediator. 
 

## Dependencies
* nuSQuIDs 

χaroν needs a modified version of nuSQuIDs which can be installed from 
  [`nuSQuIDS branch`](https://github.com/qrliu/nuSQuIDS) and please make sure python interface is activated.
* Python
* NumPy
* Scipy
* SymPy
* Astropy
* Matplotlib
* H5py


## Installation
Unzip the xsec.zip file in [`charon/xsec/`](charon/xsec) first. 
```
unzip ./charon/xsec/xsec.zip
```

data at production can be downloaded from [`charon data`](https://icecube.wisc.edu/~qliu/charon/). These data are supposed to be put into [`charon/data/`](charon/data) before installation.

χarον can be installed using the following command
```
python setup.py install
```

If you want to generate flux at production yourself, please go to [`production`](production/) and see [`README`](production/README.md) file there. 
 
## Examples
You can find examples of how to run χαρον in the 
[`examples`](examples) folder.

## Citation
If you use χarον in published work, please cite \\
preprint: arXiv:20xx.xxxxx.


For any questions, issues and comments, please contact qliu@icecube.wisc.edu.


ʕ•̫͡•ʕ•͓͡•ʔ-̫͡-ʕ•̫͡•ʔʔ-̫͡-ʔฅʕ•̫͡•ʕ•͓͡•ʔ-̫͡-ʕ•̫͡•ʔʔ-̫͡-ʔฅʕ•̫͡•ʕ•͓͡•ʔ-̫͡-ʕ•̫͡•ʔʔ-̫͡-ʔ
