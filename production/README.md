This is how to generate fluxes without EW correction implemented.
-----------------------------------------------

## Dependencies
* C++11 or higher
* [`Pythia 8`](http://home.thep.lu.se/Pythia/)


## Flux Generating

First do
```
make main
```
to compile.

To generate events, do 
```
./run.sh channel mass location process type Nevent bins lower_energy_bound binning_scale
(density mediator_mass seed) 
```

For a fast secluded run, do
```
./run_secluded.sh channel mass location Nevent bins lower_energy_bound binning_scale mediator_mass (seed)
```
It generates fluxes at densities from 0. to 155 g/cm^3. 


--channel <br/>
channel is the annihilation channel you want to generate:
Available ones are:
{'dd','uu','ss','cc','bb','tt','gg','WW','ZZ','HH','ee','mumu','tautau','nuenue','numunumu','nutaunutau'}
polarization is not included here.

--DM mass <br/>
DM mass in GeV.


--location <br/>
location of the DM annihilation. Available: "Sun", "Earth", "Halo" 


--process  <br/>
process "ann" or "decay"


--type <br/>
"-" or "secluded"


--Nevent <br/>
number of events you want to generate


--bins <br/>
number of energy bins


--lower_energy_bound <br/>
lowest energy (GeV) for the spectrum


--binning_scale <br/>
bin the spectrum in linear ("-") or log ("log") scale 
If the type is "secluded", we need to specify the density and mediator mass at the location of decay.


**Parameters for secluded**


-- density (g/cm^3) <br/>
density of the annihilation/decay location. 


-- mediator_mass (GeV) <br/>
mass of the mediator.


The last is an optional parameter 

--seed <br/>
seed for MC generation, default is 1. 

------------------------------------------------
Files are saved in ./location/ with name channel_dm mass_location_process_binscale-#.dat or ./secluded/ with name channel_dm mass_mediator mass_seed_density_bin scale-#.dat where 

0-nue

1-nue_bar

2-numu

3-numu_bar

4-nutau

5-nutau_bar 

The first column is energy is GeV. 

Generated tables are in ../charon/data with 10000000 events without EW correction  

Only fluxes at the production here. The flux is dN/dx or dN/dlogx per annihilation/decay depending on the way of binning.

For secluded DM case, fluxes at production at different locations with different matter densities along the line of sight are supposed to be generated.


ʕ•̫͡•ʕ•͓͡•ʔ-̫͡-ʕ•̫͡•ʔʔ-̫͡-ʔฅʕ•̫͡•ʕ•͓͡•ʔ-̫͡-ʕ•̫͡•ʔʔ-̫͡-ʔฅʕ•̫͡•ʕ•͓͡•ʔ-̫͡-ʕ•̫͡•ʔʔ-̫͡-ʔฅʕ•̫͡•ʕ•͓͡•ʔ-̫͡-ʕ•̫͡•ʔʔ-̫͡-ʔฅʕ•̫͡•ʕ•͓͡•ʔ-̫͡-ʕ•̫͡•ʔʔ-̫͡-ʔ 
