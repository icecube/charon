Author  : Q.R. Liu

Here is how to start running these scripts.
-----------------------------------------------
To generate events, just do 
```
./run.sh channel mass location process type Nevent bins lower_energy_bound binning_scale <br /> (density mediator_mass seed) 
```

--channel
channel is the annihilation channel you want to generate:
Now available ones are:
{'dd','uu','ss','cc','bb','tt','gg','WW','ZZ','mumu','tautau','nuenue','numunumu','nutaunutau'}


--mass

DM mass in GeV.


--location

location of the DM annihilation. Available: "Sun", "Earth" ,"Halo" 


--process

process "ann" or "decay"


--type

"-" or "secluded"


--Nevent

number of events you want to generate


--bins

number of energy bins


--lower_energy_bound

lowest energy (GeV) for the spectrum

--binning_scale

bin the spectrum in linear ("-") or log ("log") scale 
If the type is "secluded", we need to specify the density and mediator mass at the location of decay.

  -- density (g/cm^3)

  -- mediator_mass (GeV)

The last is an optional parameter 

--seed

seed for MC generation 

------------------------------------------------
Files are saved in ./location/ with name channel_mass_seed_location_process-#.dat where 

0-nue

1-nue_bar

2-numu

3-numu_bar

4-nutau

5-nutau_bar 

The first column is energy is GeV. 
"mass" here is either DM mass or mediator mass depending on the type. 

Generated tables are in ./data with name channel_mass.dat with 10000000 events in form 

|Enu|nu-e|nu-e-bar|nu-mu|nu-mu-bar|nu-tau|nu-tau-bar|

Only fluxes at the production here. The flux is dN/dx or dN/dlogx per annihilation/decay depending on the way of binning.
