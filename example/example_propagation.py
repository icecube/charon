import os,time
import numpy as np
import matplotlib.pyplot as plt 
from charon import propa
import charon.physicsconstants as PC
pc = PC.PhysicsConstants()


#info of DM
channel = 'WW'
mass    = 1000.
process = 'ann' #annihilation default

#info of neutrino spectrum binning. 
Emin         = 1.
Emax         = mass
nodes        = 100
bins         = 300

#info of oscillation parameters
theta_12=33.82
theta_13=8.6
theta_23=48.6
delta_m_12=7.39e-5
delta_m_13=2.528e-3
delta = 221.


#linear space of outcoming energies 
logscale = False #default
#include interactions
interactions = True #default

#construct the flux object.
Flux = propa.NuFlux(channel,mass,nodes,Emin=Emin,Emax=Emax,bins=bins,
                     process=process,logscale=logscale,interactions=interactions,
                     theta_12=theta_12,theta_13=theta_13,theta_23=theta_23,
                     delta_m_12=delta_m_12,delta_m_13=delta_m_13,delta=delta,pathSunModel='../charon/models/struct_b16_agss09.dat')


flux_ini_sun   = Flux.iniFlux('Sun')
flux_ini_earth = Flux.iniFlux('Earth')
flux_ini_Halo  = Flux.iniFlux('Halo')


x = Flux.iniE()/mass
plt.plot(x,flux_ini_sun['nu_mu'],linewidth=2.0,label=r'${\rm{Sun}}\;\nu_\mu$')
plt.plot(x,flux_ini_earth['nu_mu'],linewidth=2.0,label=r'${\rm{Earth}}\;\nu_\mu$')
plt.plot(x,flux_ini_Halo['nu_mu'],linewidth=2.0,label=r'${\rm{Halo}}\;\nu_\mu$')

plt.ylim(1e-4,1e3)
plt.xlim(0.0,1.0)
plt.yscale('log')
plt.xlabel(r"$x = E/m$")
plt.ylabel(r"$dN_\nu/dx\;[ann^{-1}]$")
plt.legend()
plt.savefig('/home/egenton/charon/example/initial_fluxes.png')

t0 = time.time()
flux_det_sun = Flux.Sun('detector',zenith=np.pi/6.,avg=True) 
#zenith can be replaced by latitude of detector in degree and mjd, e.g. lat_det = -90, mjd = 59062.
print ('from Sun center to Earth surface, time duration:')
print (time.time()-t0, 's')

t0 = time.time()
flux_det_Earth = Flux.Earth('detector',avg=True)
print ('from Earth center to Earth surface, time duration:')
print (time.time()-t0, 's')

t0 = time.time()
flux_det_Halo = Flux.Halo('detector',zenith=np.pi)
print ('averaged propagation, time duration:')
print (time.time()-t0, 's')