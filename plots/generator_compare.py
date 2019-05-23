#! /usr/bin/env python
#author : Q. R. Liu
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('paper.mplstyle')
from matplotlib.offsetbox import AnchoredText
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
#from hepmc import *

p = argparse.ArgumentParser(description="compare different Monte Carlo", formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--channel",default='bb',type=str,
		help='annihilation channel')
p.add_argument('--mass',default=1000,type=int,
		help='WIMP mass')
p.add_argument('--flavor',default=1,type=int,
		help='nu_e:0,nu_mu:1,nu_tau:2')
	
args = p.parse_args()

channel = args.channel
mass    = args.mass
flavor  = args.flavor


fig = plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec (4, 1,height_ratios=[4,1,1,1])
ax = fig.add_subplot(gs[0])
ax.set_yscale('log')
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)

ax.set_ylabel(r'$(1/N_{event})\:dN/dx$', fontsize=15)

ch_flavor = {1 : r'$\nu_e$', 2 : r'$\nu_\bar{e}$', 3 : r'$\nu_\mu$' , 4 : r'$\nu_\bar{\mu}$', 5 : r'$\nu_\tau$', 6 : r'$\nu_\bar{\tau}$'}
valx = np.linspace(0.,1.,201)
valx = (valx[1:]+valx[:-1])/2.

ch_wimpsim = {'dd':1,'uu':2,'ss':3,'cc':4,'bb':5,'tt':6,'gg':7,'WW':8,'ZZ':9,'mumu':10,'tautau':11,'nuenue':12,'numunumu':13,'nutaunutau':14}
wimpsim_data = np.genfromtxt("../WimpSim/data/dat-wa/wa-m1000-ch{0}-sun-sum.dat".format(ch_wimpsim[channel]))
nu = wimpsim_data[2*flavor,:]
nu_bar = wimpsim_data[2*flavor+1,:]

ax.set_title(ch_flavor[2*flavor+1]+' with '+r'$\mathrm{M}_\mathrm{DM}=\;$'+str(mass)+' GeV',fontsize=18)
ax.step(
    valx,nu, alpha=1.0, where='mid',
    linewidth=2, linestyle='-',color='steelblue',label=r'$\mathrm{WimpSim}\;$')
    
    
#old pppc
pppc_ch = {'bb':2,'tautau':3,'WW':8}

flux = []
def fit_spectra(x,a0,a1,a2,a3,a4,a5,b,c0,c1,c2):
	w = np.log10(x)
	return a0*(1+a1*w+a2*w**2+a3*w**3+a4*w**4+a5*w**5)*(1-x)**b+c0*x**c1*(1-x)**c2

if pppc_ch[channel] is not 8:
	params = np.genfromtxt("../PPPC/old/DMnuProdParameters/parameters_fit_numu_sun.dat",invalid_raise=False,dtype=None)
	for j in range(11):
		if  params[j+(pppc_ch[channel]-2)*11][0]    == mass:
			param_space =  params[j+(pppc_ch[channel]-2)*11][1:] 
	flux = map(lambda i: fit_spectra(i, *param_space), valx)
	ax.step(valx,flux,linewidth=2,color='orange',label=r'$\mathrm{Cirelli}\;\mathrm{2005}\;$')
else:
	pppc = np.genfromtxt("../PPPC/old/DMnuProdFluxes/sun_numu.dat",invalid_raise=False,dtype=None)
	x = []
	for i in range(len(pppc)):
		if pppc[i][0]== mass:
			x.append(pppc[i][1])
			flux.append(pppc[i][8])
	f = interp1d(x,flux)
	flux = f(valx)
	ax.step(valx,flux,linewidth=2,color='orange',label=r'$\mathrm{Cirelli}\;\mathrm{2005}\;$')

	

#old pppc

plot = np.genfromtxt('/data/user/qliu/DM/pythia8240/final/results/{channel}_{mass}.dat'.format(channel=channel,mass=mass))
valx = plot[:,0]/mass
valy1 = plot[:,2*flavor+1]*mass
valy2 = plot[:,2*flavor+2]*mass
ax.step(valx, valy1,alpha=1.0,where='mid',
        linewidth=2, linestyle='-', color='crimson', label=r'$\mathrm{This Work}\;$')
ax.set_xlim(0.,1.)
#ax.set_ylim(1e-4,10.0) 
ax.set_ylim(1e-6,10.0)  #bb 

for xmaj in ax.xaxis.get_majorticklocs():
    ax.axvline(x=xmaj, ls=':', color='gray', alpha=0.7, linewidth=1)
for ymaj in ax.yaxis.get_majorticklocs():
    ax.axhline(y=ymaj, ls=':', color='gray', alpha=0.7, linewidth=1)

s = {'dd':r'$d\bar{d}$','uu':r'$u\bar{u}$','ss':r'$s\bar{s}$','cc':r'$c\bar{c}$','bb':r'$b\bar{b}$','tt':r'$t\bar{t}$','gg':r'gg','WW':r'W+ W-','ZZ':r'ZZ','mumu':r'$\mu^+\mu^-$','tautau':r'$\tau^+\tau^-$','nuenue':r'$\nu_e\bar{\nu}_{e}$','numunumu':r'$\nu_\mu\bar{\nu}_{\mu}$','nutaunutau':r'$\nu_\tau\bar{\nu}_\tau$'}

at = AnchoredText(s[channel], prop=dict(size=20), frameon=True,
                  loc='lower left')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.5")
ax.add_artist(at)
ax.legend(loc='best',prop=dict(size=18),frameon=False)
#ax.legend(loc='lower right',prop=dict(size=18),frameon=False)
yticks = ax.yaxis.get_ticklabels()
yticks[0].set_visible(False)
yticks[1].set_visible(False)
ax.get_xaxis().set_ticks([])
ax.tick_params(axis='y',direction='in',which='both',left=True)

ax1 = fig.add_subplot (gs[1])
ratio = np.divide(nu,flux,out=np.zeros_like(nu),where=flux!=0)
ax1.step(valx,ratio,where='mid',linewidth=2,linestyle='-',color='olivedrab',alpha=0.7,label=r'$\mathrm{WimpSim/Cirelli\;2005}\;$')
ax1.set_ylabel(r'Ratio',fontsize=15)
ax1.legend(loc='upper right',prop=dict(size=12),frameon=False)
for ymaj in np.arange(0.,2.,0.5):
    ax1.axhline(y=ymaj, ls=':', color='gray', alpha=0.7, linewidth=1)
ax1.axhline(y=1., ls=':', color='black', alpha=1, linewidth=1.5)
ax1.set_xlim(0.,1.)
ax1.set_ylim(0.,2.)
ax1.set_yticks([0.5,1,1.5])
ax1.tick_params(axis='both',direction='in',which='both',bottom=True,left=True)

ax2 = fig.add_subplot (gs[2])
ratio = np.divide(valy1,flux,out=np.zeros_like(valy1),where=flux!=0)
ax2.step(valx,ratio,where='mid',linewidth=2,linestyle='-',color='chocolate',alpha=0.7,label=r'$\mathrm{This\;Work/Cirelli\;2005}\;$')
ax2.set_ylabel(r'Ratio',fontsize=15)
ax2.legend(loc='upper right',prop=dict(size=12),frameon=False)
for ymaj in np.arange(0.,2.,0.5):
    ax2.axhline(y=ymaj, ls=':', color='gray', alpha=0.7, linewidth=1)
ax2.axhline(y=1., ls=':', color='black', alpha=1, linewidth=1.5)
ax2.set_xlim(0.,1.)
ax2.set_ylim(0.,2.)
ax2.set_yticks([0.5,1,1.5])
ax2.tick_params(axis='both',direction='in',which='both',bottom=True,left=True)

ax3 = fig.add_subplot (gs[3])
ratio = np.divide(valy1,nu,out=np.zeros_like(valy1),where=nu!=0)
ax3.step(valx,ratio,where='mid',linewidth=2,linestyle='-',color='magenta',alpha=0.7,label=r'$\mathrm{This\;Work/WimpSim}\;$')
ax3.set_ylabel(r'Ratio',fontsize=15)
#ax3.legend(loc='best',prop=dict(size=12),frameon=False)
ax3.legend(loc='upper right',prop=dict(size=12),frameon=False)
for ymaj in np.arange(0.,2.,0.5):
    ax3.axhline(y=ymaj, ls=':', color='gray', alpha=0.7, linewidth=1)
ax3.axhline(y=1., ls=':', color='black', alpha=1, linewidth=1.5)
ax3.set_xlim(0.,1.)
ax3.set_ylim(0.,2.)
ax3.set_yticks([0.5,1,1.5])
ax3.tick_params(axis='both',direction='in',which='both',bottom=True,left=True)


plt.xlabel(r'$\mathrm{x}=\mathrm{E}/\mathrm{M}_{\mathrm{DM}}$', fontsize=18)
 
fig.savefig('{channel}_{mass}_{flavor}.pdf'.format(channel=channel, mass=mass, flavor=flavor), bbox_inches='tight', dpi=150)

