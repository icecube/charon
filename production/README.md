Date    : May 24 2019
Author  : Q.R. Liu

Here is how to start running these scripts.
-----------------------------------------------
To generate events, just do 
```
./run.sh channel mass bins seed Nevent location
```
--channel

channel is the annihilation channel you want to generate:
Now available ones are:
{'dd','uu','ss','cc','bb','tt','gg','WW','ZZ','mumu','tautau','nuenue','numunumu','nutaunutau'}

--mass

DM mass in GeV.

--bins

Number of energy bins

--Nevent

number of events you want to generate.

--location

location of the DM annihilation. Available: "Sun", "Earth" ,"Galactic" 

--process

process "ann" or "decay"

--type

"-" or "secluded"

if the type is "secluded", we need to specify the density and mediator mass at the location of decay.
  
  -- density

  -- mediator mass

The last is an optional parameter 

--seed

seed for MC generation 


------------------------------------------------
Files are saved in ./location/ with name channel_mass_seed_location-#.dat where
 
0-nue

1-nue_bar

2-numu

3-numu_bar

4-nutau

5-nutau_bar 

and the first column is energy is GeV.

Generated tables are in ./data with name channel_mass.dat with 10000000 events in form 
|Enu|nu-e|nu-e-bar|nu-mu|nu-mu-bar|nu-tau|nu-tau-bar|

Only fluxes at the production here. The flux is dN/dE per annihilation.
