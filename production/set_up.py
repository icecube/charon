#initialize settings of the generator
#Author  :  Q.R. Liu

import numpy as np
import argparse
import os 
import sympy as sym
from sympy import Symbol
from sympy.solvers import solve

p = argparse.ArgumentParser(description="parameters needed to make .cmnd files", formatter_class=argparse.RawTextHelpFormatter)

p.add_argument("--channel",default='bb',type=str,
				help="channel of the DM annihilation")
p.add_argument("--process",default='ann',type=str,
				help="annihilation or decay")
p.add_argument("--mass",default="1000.",type=str,
				help="DM mass")
p.add_argument("--mass_phi",default="1000.",type=str,
				help="mediator mass")
p.add_argument("--bins",default=200,type=int,
				help="binning")
p.add_argument("--Nevent",default=1000,type=int,
				help="number of events")
p.add_argument("--location",default='Sun',type=str,
				help="Location of DM annihilation")
p.add_argument("--type",default='-',type=str,
				help="Is it secluded DM ? '-' or 'secluded'")
args = p.parse_args()

channel      = args.channel
str_mass	 = args.mass
bins	     = args.bins
Nevent	     = args.Nevent
location     = args.location
process      = args.process
process_type = args.type
str_mphi     = args.mass_phi
me = 0.000511

mass = float(str_mass)
mphi = float(str_mphi)

if   process == 'ann':
	Ecm     = 2*mass
elif process == 'decay':
	Ecm     = mass

	
code    = {'dd':1,'uu':2,'ss':3,'cc':4,'bb':5,'tt':6,'gg':21,'ZZ':23,'WW':24,'HH':25,'ee':11,'nuenue':12,'mumu':13,'numunumu':14,'tautau':15,'nutaunutau':16}
nu_number  = {'nuenue':0,'numunumu':2,'nutaunutau':4}

if os.path.exists('./output')   ==False:
	os.mkdir('./output')
if os.path.exists('./cmnd')     ==False:
	os.mkdir('./cmnd')
if os.path.exists('./'+location)==False:
	os.mkdir('./'+location)


if process_type == 'secluded':
	if os.path.exists('./secluded')==False:
		os.mkdir('./secluded')
	Ecm     = mass

if channel in ('nuenue','numunumu','nutaunutau'):
	binning    = np.linspace(0.,mass,bins+1)
	bin_center = (binning[1:]+binning[:-1])/2.
	for i in range(6):
		flux_path = './{}/{}_{}_{}_{}-{}.dat'.format(location,channel,mass,location,process,i)	
		f = open(flux_path,'w')
		if i == nu_number[channel] or i==nu_number[channel]+1 : 
			diff = np.diff(binning)	
			if process == 'ann':
				for j in range(bins):
					if j < bins-1:  
						f.write('{:.4E}'.format(bin_center[j])+'\t'+'0.0000e+00'+'\n')
					else:
						non_zero = 1./diff[-1]	
						f.write('{:.4E}'.format(bin_center[j])+'\t'+'{}'.format(non_zero)+'\n')
			elif process == 'decay':
				if (bins % 2) == 0:
					mid      = bin_center[bins/2-1] 
					non_zero = diff[bins/2-1] 
				else:
					mid      = bin_center[(bins-1)/2] 
					non_zero = diff[(bins-1)/2] 

				for j in range(bins):
					if  bin_center[j] == mid:
						f.write('{:.4E}'.format(bin_center[j])+'\t'+'{}'.format(non_zero)+'\n')
					else:
						f.write('{:.4E}'.format(bin_center[j])+'\t'+'0.0000e+00'+'\n')


		else:
			for j in range(bins):
				f.write('{:.4E}'.format(bin_center[j])+'\t'+'0.0000e+00'+'\n')

	        
		f.close()


else:
	if process_type == 'secluded':
		f_name  ='./cmnd/'+channel+'_'+str_mass+'_'+str_mphi+'.cmnd'
	else:
		f_name  ='./cmnd/'+channel+'_'+str_mass+'_'+process+'.cmnd'
	
	f = open(f_name,'w')
	f.write('Main:numberOfEvents = '+str(Nevent)+'\n')
	f.write('Main:timesAllowErrors = 5'+'\n') 
	f.write('Init:showChangedSettings = on'+'\n')     
	f.write('Init:showChangedParticleData = on'+'\n') 
	f.write('Next:numberCount = 100'+'\n')            
	f.write('Next:numberShowInfo = 1'+'\n')           
	f.write('Next:numberShowProcess = 1'+'\n')        
	f.write('Next:numberShowEvent = 1'+'\n') 
	f.write('Beams:idA = -11'+'\n')   
	f.write('Beams:idB = 11'+'\n')    
	f.write('PDF:lepton = off'+'\n')
	if process_type == 'secluded':
		x = Symbol('x')
		y = Symbol('y')
		px= Symbol('px')
		py= Symbol('py')
		solution = solve([x+y-mass,px+py-np.sqrt(mass**2-mphi**2),x**2-px**2-me**2,y**2-py**2-me**2],dict=True)
		eA = float(solution[0][y])
		eB = float(solution[0][x])
		f.write('Beams:frameType = 2'+'\n')
		f.write('Beams:eA = '+str(eA)+'\n')
		f.write('Beams:eB = '+str(eB)+'\n')
		f.write('999999:all = GeneralResonance void 1 0 0 '+str(mphi)+' 1. 0. 0. 0.'+'\n')
	
	else:
		f.write('Beams:eCM = '+str(Ecm)+'\n')
		f.write('999999:all = GeneralResonance void 1 0 0 '+str(Ecm)+' 1. 0. 0. 0.'+'\n')
	if code[channel] in [21,22,23,25]:
		f.write('999999:addChannel = 1 1.0 101 '+str(code[channel])+' '+str(code[channel])+'\n') 
	else:
		f.write('999999:addChannel = 1 1.0 101 '+str(code[channel])+' '+str(-code[channel])+'\n')    
	
	f.close()
