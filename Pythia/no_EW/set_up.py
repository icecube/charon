#initialize settings of the generator
#Author  :  Q.R. Liu

import numpy as np
import argparse

p = argparse.ArgumentParser(description="parameters needed to make .cmnd files", formatter_class=argparse.RawTextHelpFormatter)

p.add_argument("--channel",default='bb',type=str,
				help="channel of the DM annihilation")
p.add_argument("--mass",default=1000,type=int,
				help="DM mass")
p.add_argument("--Nevent",default=1000,type=int,
				help="number of events")
args = p.parse_args()

channel = args.channel
mass    = args.mass
Ecm     = float(2*mass)
Nevent  = args.Nevent

code    = {'dd':1,'uu':2,'ss':3,'cc':4,'bb':5,'tt':6,'gg':21,'WW':24,'ZZ':9,'mumu':13,'tautau':15,'nuenue':12,'numunumu':14,'nutaunutau':16}

f = open(channel+'_'+str(mass)+'.cmnd','w')
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
f.write('Beams:eCM = '+str(Ecm)+'\n')

f.write('999999:all = GeneralResonance void 1 0 0 '+str(Ecm)+' 1. 0. 0. 0.'+'\n')
f.write('999999:addChannel = 1 1.0 101 '+str(code[channel])+' '+str(-code[channel])+'\n')    

f.close()
