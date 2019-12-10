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
{'dd':1,'uu':2,'ss':3,'cc':4,'bb':5,'tt':6,'gg':7,'WW':8,'ZZ':9,'mumu':10,'tautau':11,'nuenue':12,'numunumu':13,'nutaunutau':14}
numbers correpond to the MC particle numbering scheme

--mass
DM mass in GeV.

--bins
Number of energy bins

--seed
seed for MC generation 

--Nevent
number of events you want to generate.

--location
location of the DM annihilation. Available: "Sun","Earth" 

------------------------------------------------

Files generated are in ./data with name channel_mass.dat with 10000000 events in form 
|Enu|nu-e|nu-e-bar|nu-mu|nu-mu-bar|nu-tau|nu-tau-bar|

Only fluxes at the production here. The flus is dN/dE per annihilation.
