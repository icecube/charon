#author : Q.R. Liu

import os
import numpy as np
import scipy as sp
import scipy.special as spe
import scipy.interpolate as interpolate
from sympy.solvers import solve
from sympy import Symbol
# my modules
import neutrinocommon
import neutrinocommon.neu.xsections as xs
import neutrinocommon.physconst.physicsconstants as PC
import neutrinocommon.neu.neuosc as no
import neutrinocommon.astro.body as bd
import neutrinocommon.exp.icecube as ice
import neutrinocommon.tools.generaltools as gt
import nuSQUIDSpy as nsq
import nuSQUIDSTools
import DM

pc = PC.PhysicsConstants()

datDMFluxSweden = '/data/user/qliu/DM/dat-wa/'
#wimpsim channels

ch_wimpsim = {'dd':1,'uu':2,'ss':3,'cc':4,'bb':5,'tt':5,'gg':7,'ww':8,'ZZ':9,'mumu':10,'tautau':11,'nuenue':12,'numunumu':13,'nutaunutau':14}

neu_wimpsim = {1 : '$\nu_e$', 2 : '$\nu_\bar{e}$', 3 : '\nu_\mu' , 4 : '$\nu_\bar{\mu}$', 5 : '$\nu_\tau$', 6 : '$\nu_\bar{\tau}$'}

#wimpsim
def DMSweFlux(Enu,neuflavor,ch,DMm,location = 'SunCtr'):
    """ Gets dN/dz(z) using de swedish data set.
    @type  Enu          :      float
    @param Enu          :      neutrino energy [GeV]
    @type  neuflavor    :      integer
    @param neuflavor    :      flavor of the neutrino (0 : nu_e, 1 : anu_e, 2 : nu_mu 0, 3 : anu_mu, 4 : nu_tau, 5 : anu_tau)
    @type  ch           :      string
    @param ch           :      annihilation channel
    @type  DMm          :      float
    @param DMm          :      dark matter mass [GeV]
    @type location      :      string
    @param location	:      SunCtr line 1-6, SunSrfc 7-12,Sun1Au 13-18, SunSrfc2nd 19-24, Sun1Au2nd 25-30 
    @rtype              :      array
    @return             :      z = E/DMm, dN/dz(z) arrays at production (??)
    """
    if PC.act_channel == ch and PC.act_DM_mass == DMm and PC.act_neuflavor == neuflavor and PC.flag_inter:
        if Enu/DMm < 0.0025:
            return 0.0
        elif Enu/DMm <= 0.9975:    
            return PC.act_inter(Enu/DMm)
        elif Enu/DMm > 0.9975:
            return 0.0
        else :
            print "Interpolation error."
            quit()
    else:
        #print "reloading DM initial flux"
        DMmstring = format(DMm,'.0f')
        filename = "wa-m"+str(DMmstring)+"-ch"+str(ch_wimpsim[ch])+"-sun-sum.dat"
        file = open(datDMFluxSweden + filename,'r')
        z=np.linspace(0.0,1.,200,endpoint=False)
        #h,dat = gt.hreadfilev2(file)
	if location == 'SunCtr':
		line = neuflavor
	elif location == 'SunSrfc':
		line = neuflavor+6
	elif location == 'Sun1AU':
		line = neuflavor+12
	elif location == 'SunSrfc2nd':
		line = neuflavor+18
	elif location == 'Sun1AU2nd':
		line = neuflavor+24
	else:
		print "No such annihilation location"
		quit()
	#print location
	dn_dz =  np.genfromtxt(file)[line,:]
	#xneuflavor=neuflavor
        #if AtEarth:
        #    xneuflavor+=12
        #dn_dz = dat[xneuflavor]
        #print dat
        
        #PC.act_channel, PC.act_DM_mass, PC.act_neuflavor,PC.flag_inter = ch,DMm,neuflavor,'SunCtr'
        
        if Enu/DMm < z[0]:
            return 0.0
        elif Enu/DMm <= z[-1]:    
            inter = sp.interpolate.interp1d(z,dn_dz)
            #inter = sp.interpolate.UnivariateSpline(z,dn_dz)
            PC.act_inter = inter
            PC.flag_inter = True
            return inter(Enu/DMm)
        elif Enu/DMm > z[-1]:
            return 0.0
        else :
            print "Interpolation Error."
            quit()
	file.close()
## Creating a DM distribution ##

class DM_distribution():
    
    def __init__(self,ch,DMm,flavor):
        """ Initializes DM distribution for a given channel and "flavor".
    
        @type  ch            :      string
        @param ch            :      annihilation channel
        @type  DMm           :      float
        @param DMm           :      dark matter mass [GeV]
        @type  flavor        :      integer
        @param flavor        :      flavor of the neutrino (0 : nu_e, 1 : anu_e, 2 : nu_mu 0, 3 : anu_mu, 4 : nu_tau, 5 : anu_tau)
        """
        self.min    = 0.0
        self.max    = DMm
        self.DMm    = DMm
        self.ch     = ch
        self.flavor = flavor
        
    def PDF(self,Enu,location):
        """ Calculates dN/dz where z = E_nu/X_mass using the swedish data.
    
        @type  Enu          :      float
        @param Enu          :      neutrino energy [GeV]
    
        @rtype              :      float
        @return             :      dn/dz [dimensionless]
        """
        return DMSweFlux(Enu,self.flavor,self.ch,self.DMm,location)#/self.DMm

def DMSweFluxEarth(Enu,neuflavor,ch,DMm,DMsig,unit,param):
    """ Gets dN/dz(z) using de swedish data set.
    @type  Enu          :      float
    @param Enu          :      neutrino energy [GeV]
    @type  neuflavor    :      integer
    @param neuflavor    :      flavor of the neutrino (0 : nu_e, 1 : anu_e, 2 : nu_mu , 3 : anu_mu, 4 : nu_tau, 5 : anu_tau)
    @type  ch           :      string
    @param ch           :      annihilation channel
    @type  DMm          :      float
    @param DMm          :      dark matter mass [GeV]
    @type  unit         :      string
    @param unit 	:      'WimpAnn' or flux
    @rtype              :      array
    @return             :      z = E/DMm, dN/dz(z) arrays at production (??)
    """
    if PC.act_channel == ch and PC.act_DM_mass == DMm and PC.act_neuflavor == neuflavor and PC.flag_inter:
        if Enu/DMm < 0.0025:
            return 0.0
        elif Enu/DMm <= 0.9975:    
            return PC.act_inter(Enu/DMm)
        elif Enu/DMm > 0.9975:
            return 0.0
        else :
            print "Interpolation error."
            quit()
    else:
        #print "reloading DM initial flux"
        DMmstring = format(DMm,'.0f')
        filename = "wa-m"+str(DMmstring)+"-ch"+str(ch_wimpsim[ch])+"-sun-sum.dat"
        file = open(datDMFluxSweden + filename,'r')
        bins = np.linspace(0.0,1.,201,endpoint=True)
	z    = (bins[1:]+bins[:-1])/2.
        #h,dat = gt.hreadfilev2(file)
	#line = neuflavor+12
	line = neuflavor+6
	dn_dz =  np.genfromtxt(file)[line,:]
	#print dn_dz
	
	if unit == 'WimpAnn':
            if Enu/DMm < z[0]:
                return 0.0
            elif Enu/DMm <= z[-1]:    
                inter = sp.interpolate.interp1d(z,dn_dz)
                PC.act_inter = inter
                PC.flag_inter = True
                return inter(Enu/DMm)
            elif Enu/DMm > z[-1]:
                return 0.0
            else :
                print "Interpolation Error."
                quit()
	    file.close()
	elif unit == 'flux':
	    z = z*DMm
	    DM_annihilation_rate_Sun = float(np.sum(DM.DMSunAnnihilationRate(DMm*param.GeV,DMsig,param)))*param.sec # ann/s	
	    #factor = DM_annihilation_rate_Sun/(4.0*np.pi*(param.AU/param.cm)**2*DMm)
	    factor = DM_annihilation_rate_Sun/(4.0*np.pi*(param.SUNRADIUS*param.km/param.cm)**2*DMm)
	    dn_dz = factor * dn_dz
            if Enu < z[0]:
                return 0.0
            elif Enu <= z[-1]:    
                inter = sp.interpolate.interp1d(z,dn_dz)
            #inter = sp.interpolate.UnivariateSpline(z,dn_dz)
                PC.act_inter = inter
                PC.flag_inter = True
		#print inter(Enu)
                return inter(Enu)
            elif Enu > z[-1]:
                return 0.0
            else :
                print "Interpolation Error."
                quit()
  	    file.close()
	else:
	   print "No such data"
	   quit()
        

##END SWEDISH WAY
#pythia8
def pythiaflux(Enu,neuflavor,ch,DMm):
    flavor = {0 : 'e', 1 : 'ebar', 2 : 'mu', 3 : 'mubar', 4 : 't', 5 :'tbar'}

    if PC.act_channel == ch and PC.act_DM_mass == DMm and PC.act_neuflavor == neuflavor and PC.flag_inter:
        if Enu/DMm < 0.0025:
            return 0.0
        elif Enu/DMm <= 0.9975:    
            return PC.act_inter(Enu/DMm)
        elif Enu/DMm > 0.9975:
            return 0.0
        else :
            print "Interpolation error."
            quit()
    else:
    	filedir = '/data/user/qliu/DM/pythia_sim/'
  	DMmstring = format(DMm,'.0f')
   	#filename = "nu_"+str(DMmstring)+'_'+ch+'-'+flavor[neuflavor]+".dat"
	#E = np.genfromtxt(filedir+filename)[:,0]
	#dNdx = np.genfromtxt(filedir+filename)[:,1]*DMm
	if ch=='ww':
		filename = 'WW'+'_'+str(DMmstring)+".dat"
	else:
		filename = ch+'_'+str(DMmstring)+".dat"
	E = np.genfromtxt(filedir+filename)[:,0]
	dNdx = np.genfromtxt(filedir+filename)[:,1+neuflavor]*DMm
        
	if Enu < E[0]:
            return 0.0
        elif Enu <= E[-1]:    
            inter = sp.interpolate.interp1d(E,dNdx)
            PC.act_inter = inter
            PC.flag_inter = True
            return inter(Enu)
        elif Enu > E[-1]:
            return 0.0
        else :
            print "Interpolation Error."
            quit()


def herwigflux(Enu,neuflavor,ch,DMm):
    flavor = {0 : 12, 1 : -12, 2 :14, 3 : -14, 4 :16, 5 :-16}
    flavor_name = {0 : 'nu_e', 1 : 'nu_ebar', 2 :'nu_mu', 3 :'nu_mubar', 4 :'nu_tau', 5 :'nu_taubar'}
    filedir = '/data/user/qliu/DM/herwig/'

    if PC.act_channel == ch and PC.act_DM_mass == DMm and PC.act_neuflavor == neuflavor and PC.flag_inter:
        if Enu/DMm < 0.0025:
            return 0.0
        elif Enu/DMm <= 0.9975:    
            return PC.act_inter(Enu/DMm)
        elif Enu/DMm > 0.9975:
            return 0.0
        else :
            print "Interpolation error."
            quit()
    else:
  	DMmstring = format(DMm,'.0f')
   	#data = np.load(filedir+ch+'_'+str(DMmstring)+"_"+flavor_name[neuflavor]+".npy")
   	data = np.load(filedir+ch+'-'+str(DMmstring)+"_"+flavor_name[neuflavor]+".npy")
	E = np.linspace(0.,DMm,201)
	#E = np.linspace(0.,DMm,11)
        
	#event nunmber
	#ww
	evtnum = 127566
	#evtnum = 10000
    	bin_center = (E[1:]+E[:-1])/2.
   	bin_width  = np.diff(E)[0]
    	weight     = 1./(evtnum*bin_width)
    	nu_hist, _ = np.histogram(data,bins = E,weights=[weight]*len(data))
    	dNdx       = nu_hist*DMm
	#print 'Enu'
	#print Enu
	#print 'bin_center'
	#print bin_center
	if Enu < bin_center[0]:
            return 0.0
        elif Enu <= bin_center[-1]:    
            inter = sp.interpolate.interp1d(bin_center,dNdx)
            PC.act_inter = inter
            PC.flag_inter = True
            return inter(Enu)
        elif Enu > bin_center[-1]:
            return 0.0
        else :
            print "Interpolation Error."
            quit()


def SunZenith(MJD,l_det):
	#Sun Zenith in radian
	JD= MJD+2400000.5
	n = JD-2451545.
	L = 280.460+0.9856474*n
	g = 357.528+0.9856003*n
	Lambda = L+1.915*np.sin(np.deg2rad(g))+0.020*np.sin(np.deg2rad(2*g))
	number = Lambda//360
	Lambda = Lambda - number*360.
	epsilon = 23.439-0.0000004*n
	delta = np.arcsin(np.sin(np.deg2rad(epsilon))*np.sin(np.deg2rad(Lambda)))
	return abs(delta-np.deg2rad(l_det))

def Distance(theta,param):
	#distance of total vacuum, atmosphere+earth in km
	#input zenith in radian
	AU       = param.AU/param.km
	r_earth  = param.EARTHRADIUS
	x        = Symbol('x')
	solution = solve((x*x+r_earth*r_earth-AU*AU)/(2*x*r_earth)-np.cos(np.pi-theta),x)
	d_tot     = [f for f in solution if f > 0][0]
	#if np.pi/2.< theta <= np.pi:
	#	d_earth  = 2*np.cos(np.pi-theta)*param.EARTHRADIUS 
	#else:
	#	d_earth  = 0.
	d_earthatm    = nsq.EarthAtm.Track(theta)
	d_earthatm    = d_earthatm.GetFinalX()/param.km
	d_vacuum = d_tot-d_earthatm
	print d_tot, d_vacuum, d_earthatm
	return np.array([d_tot,d_vacuum, d_earthatm])


def NuFluxSurface_multi_both(proflux,Enu_min,Enu_max,nodes,ch,DMm, DMsig,param,theta_12=33.82,theta_23=49.7,theta_13=8.61,delta_m_12=7.39e-5,delta_m_13=2.525e-3,delta=0.,logscale=False,interactions=False,location = 'Earth',time=57754.,angle=None,latitude=-90.):
	''' calculate neutrino flux at sun surface (multiple energy mode with interactions)
	@type  Enu	:	float
	@param Enu	:	GeV	
	@type  DMm	:	float
	@param DMm	:	GeV
	@param theta    :       zenith
	@return 	:	flux GeV^-1sr^-1s^-1cm^-2
	'''
	#annihilation
	
	#DM_annihilation_rate_Sun = float(np.sum(DM.DMSunAnnihilationRate(DMm*param.GeV,DMsig,param)))*param.sec
	DM_annihilation_rate_Sun = 1.
	E_nodes = nodes
	Enu_min = Enu_min*param.GeV
	Enu_max = Enu_max*param.GeV
	#e_range = np.linspace(Enu_min*pc.GeV,Enu_max*pc.GeV,100)
	if logscale is True:
		e_vector = np.logspace(np.log10(Enu_min),np.log10(Enu_max),E_nodes)
	else:
		e_vector = np.linspace(Enu_min,Enu_max,E_nodes)
	

	
	if proflux == 'Pythia':
	    production = pythiaflux
	elif proflux == 'Herwig':
            production = herwigflux
	elif proflux == 'WimpSim':
	    production = DMSweFlux
   	else:	
            print 'No Such Production'

	#factor = 1  
	flux = {}
	#xsec = nsq.NeutrinoDISCrossSectionsFromTables('./xsec/nusigma_')
	#xsec = nsq.NeutrinoDISCrossSectionsFromTables('/data/user/qliu/DM/GOLEMTools/sources/nuSQuIDS/data/xsections/nusigma_')
	#nuSQ = nsq.nuSQUIDS(e_vector,3,nsq.NeutrinoType.both,interactions,xsec)
	nuSQ = nsq.nuSQUIDS(e_vector,3,nsq.NeutrinoType.both,interactions)
	energy = nuSQ.GetERange()
	for i in range(3):
		flux[str(i)+'_nu'] = np.array(map(lambda E_nu: production(E_nu/param.GeV,i*2,ch,DMm),energy))
		flux[str(i)+'_nubar'] = np.array(map(lambda E_nu: production(E_nu/param.GeV,i*2+1,ch,DMm),energy))
	print 'E_range'
	print energy 
	print 'flux'
	print flux
	#nuSQ.Set_Body(nsq.ConstantDensity(13.0,0.5))
	nuSQ.Set_Body(nsq.Sun())
	nuSQ.Set_Track(nsq.Sun.Track(param.SUNRADIUS*param.km))
	#nuSQ.Set_Track(nsq.ConstantDensity.Track(param.SUNRADIUS*param.km))
	#nuSQ.Set_MixingParametersToDefault()
	nuSQ.Set_MixingAngle(0,1,np.deg2rad(theta_12))
	nuSQ.Set_MixingAngle(0,2,np.deg2rad(theta_13))
	nuSQ.Set_MixingAngle(1,2,np.deg2rad(theta_23))
	#nuSQ.Set_MixingAngle(0,3,np.deg2rad(theta_23))
	nuSQ.Set_SquareMassDifference(1,delta_m_12)
	nuSQ.Set_SquareMassDifference(2,delta_m_13)
	nuSQ.Set_abs_error(1.e-5)
	nuSQ.Set_rel_error(1.e-5)
	nuSQ.Set_CPPhase(0,2,delta)
	nuSQ.Set_ProgressBar(True)
	initial_flux = np.zeros((E_nodes,2,3))
	for j in range(len(flux['0_nu'])):
		for k in range(3):
			initial_flux[j][0][k] = flux[str(k)+'_nu'][j]
			initial_flux[j][1][k] = flux[str(k)+'_nubar'][j]
	print 'initial_flux'
	print initial_flux
	nuSQ.Set_initial_state(initial_flux,nsq.Basis.flavor)
	nuSQ.Set_TauRegeneration(True)
	nuSQ.EvolveState()
	
	e_range = np.linspace(Enu_min,Enu_max,nodes)
	flux_surface = np.zeros(len(e_range),dtype = [('nu_e','float'),('nu_mu','float'),('nu_tau','float'),('nu_e_bar','float'),('nu_mu_bar','float'),('nu_tau_bar','float'),('zenith','float')]) 
	

	if location == 'Sunsfc':
		#factor = DM_annihilation_rate_Sun/(4.0*np.pi*(param.SUNRADIUS*param.km/param.cm)**2*DMm)
		factor = 1. 
		flux_surface['nu_e'] = factor*np.array([nuSQ.EvalFlavor(0,e,0) for e in   e_range])
		flux_surface['nu_mu'] = factor*np.array([nuSQ.EvalFlavor(1,e,0) for e in  e_range])
		flux_surface['nu_tau'] = factor*np.array([nuSQ.EvalFlavor(2,e,0) for e in e_range])
		flux_surface['nu_e_bar'] = factor*np.array([nuSQ.EvalFlavor(0,e,1) for e in e_range])
		flux_surface['nu_mu_bar'] = factor*np.array([nuSQ.EvalFlavor(1,e,1) for e in e_range])
		flux_surface['nu_tau_bar'] = factor*np.array([nuSQ.EvalFlavor(2,e,1) for e in e_range])
	elif location == 'Earth':
		if angle == None:
			zenith = SunZenith(time,latitude)
		else:
			print (angle)
			zenith = np.deg2rad(angle)
		d_tot, d_vacuum, d_earthatm = Distance(zenith,param)
		
		#factor = DM_annihilation_rate_Sun/(4.0*np.pi*(param.AU/param.cm)**2*DMm)
		factor = DM_annihilation_rate_Sun/(4.0*np.pi*(d_tot*param.km/param.cm)**2*DMm)

	
		composition_new =np.array([[[nuSQ.EvalFlavor(0,e,0),nuSQ.EvalFlavor(1,e,0),nuSQ.EvalFlavor(2,e,0)],[nuSQ.EvalFlavor(0,e,1),nuSQ.EvalFlavor(1,e,1),nuSQ.EvalFlavor(2,e,1)]] for e in e_range])
	#print composition_new
		
		nuSQ.Set_Body(nsq.Vacuum())
		
		nuSQ.Set_Track(nsq.Vacuum.Track(float(d_vacuum)*param.km))
		nuSQ.Set_ProgressBar(True)
	
		nuSQ.Set_initial_state(composition_new,nsq.Basis.flavor)
		nuSQ.Set_TauRegeneration(True)
		nuSQ.EvolveState()
		
		composition_new =np.array([[[nuSQ.EvalFlavor(0,e,0),nuSQ.EvalFlavor(1,e,0),nuSQ.EvalFlavor(2,e,0)],[nuSQ.EvalFlavor(0,e,1),nuSQ.EvalFlavor(1,e,1),nuSQ.EvalFlavor(2,e,1)]] for e in e_range])
		
		nuSQ.Set_Body(nsq.EarthAtm())
		nuSQ.Set_Track(nsq.EarthAtm.Track(zenith))
		nuSQ.Set_ProgressBar(True)
	
		nuSQ.Set_initial_state(composition_new,nsq.Basis.flavor)
		nuSQ.Set_TauRegeneration(True)
		nuSQ.EvolveState()
		flux_surface['nu_e']       = factor*np.array([nuSQ.EvalFlavor(0,e,0) for e in   e_range])
		flux_surface['nu_mu']      = factor*np.array([nuSQ.EvalFlavor(1,e,0) for e in  e_range])
		flux_surface['nu_tau']     = factor*np.array([nuSQ.EvalFlavor(2,e,0) for e in e_range])
		flux_surface['nu_e_bar']   = factor*np.array([nuSQ.EvalFlavor(0,e,1) for e in e_range])
		flux_surface['nu_mu_bar']  = factor*np.array([nuSQ.EvalFlavor(1,e,1) for e in e_range])
		flux_surface['nu_tau_bar'] = factor*np.array([nuSQ.EvalFlavor(2,e,1) for e in e_range])
		flux_surface['zenith']	   = np.array([zenith]*len(e_range))	
	
	return flux_surface
