"""
Author  : C.A. Arguelles
Date    : 10/MAY/2011

Implements routines to calculate dark matter annihilation rates.

Log :
- Modified on 30/MAY/12 : added references.
"""

import numpy as np
import scipy as sp
import scipy.special as spe
import scipy.interpolate as interpolate
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
# my modules
import neutrinocommon
import neutrinocommon.neu.xsections as xs
import neutrinocommon.physconst.physicsconstants as PC
import neutrinocommon.neu.neuosc as no
import neutrinocommon.astro.body as bd
import neutrinocommon.exp.icecube as ice
import neutrinocommon.tools.generaltools as gt
#import neutrinocommon.shareobjects.optneuosc as ono

pc = PC.PhysicsConstants()

filepath = neutrinocommon.astro.__path__
datDMFlux = filepath[0] + "/data/dm/DMnuProdParameters/"
#datDMFluxSweden = filepath[0] + "/data/dm/wimpsim/"
datDMFluxSweden = '/home/qliu/DM/dat-wa'

numuSun = "parameters_fit_numu_sun.dat"
nutauSun = "parameters_fit_nutau_sun.dat"

ch = {'bb' : 0, 'tautau' : 1, 'cc': 2, 'ss' : 3, 'gg' : 4}
# wimpsim names-tags
# http://copsosx03.physto.se/wimpsim/channels.html

#ch_wimpsim  = {'cc' : 1 , 'bb' : 2, 'tt' : 3, 'tautau' : 4, 'WW' : 5, 'ZZ' : 6, 'gg' : 12, 'uu' : 14, 'ss' : 15, 'nuenue' : 16, 'numunumu' : 17, 'nutaunutau' : 18, 'dd' : 19}
ch_wimpsim = {'dd':1,'uu':2,'ss':3,'cc':4,'bb':5,'tt':5,'gg':7,'WW':8,'ZZ':9,'mumu':10,'tautau':11,'nuenue':12,'numunumu':13,'nutaunutau':14}

neu_wimpsim = {1 : '$\nu_e$', 2 : '$\nu_\bar{e}$', 3 : '\nu_\mu' , 4 : '$\nu_\bar{\mu}$', 5 : '$\nu_\tau$', 6 : '$\nu_\bar{\tau}$'}
DM = {10 : 0, 30 : 1, 50: 2, 70 : 3, 90 : 4, 100 : 5, 300 : 6, 500 : 7, 700 : 8, 900 : 9, 1000 : 10}

def ChannelLabel(ch):
    """ Returns label of a given annihilation channel.
    @type  ch    :      string
    @param ch    :      annihilation channel
    
    @rtype       :      string
    @return      :      returns label for the annihilation channel
    """
    if ch == 'bb' :
        return 'b\\bar{b}'
    elif ch == 'tautau' :
        return '\\tau\\bar{\\tau}'
    elif ch == 'cc' :
        return 'c\\bar{c}'
    elif ch == 'ss' :
        return 's\\bar{s}'
    elif ch == 'gg' :
        return 'gg'
    else :
        print 'Unavailable channel label.'
        quit()

##BEGIN ITALIAN WAY

def DMParameters(neuflavor):
    """ Gets DM parameters for a given neutrino flavor.
    @type  neuflavor    :      integer
    @param neuflavor    :      flavor of the neutrino (0 : e, 1 : mu, 2 : tau)
    
    @rtype              :      array
    @return             :      returns the parametrization of the neutrino flux at production.
    """
    if neuflavor == 0 or neuflavor == 1 :
        file = open(datDMFlux + numuSun,'r')
    elif neuflavor == 2 :
        file = open(datDMFlux + nutauSun,'r')
    else:
        print "Unsupported neutrino flavor in DM parameters."
        quit()
    h,p = gt.hreadfilev3(file)
    return p

def DMFlux(Enu,DMm,c,p):
    """ Calculates DMm in channel 'c' according to reference
    
    Ref. ec. 7 arXiv : hep-ph/0506298
    
    @type   Enu     :  float
    @param  Enu     :  neutrino energy [eV]
    @type   DMm     :  float
    @param  DMm     :  dark matter mass [eV]
    @type   c       :  string
    @param  c       :  channel 
    @type   p       :  list
    @param  p       :  dark matter parameters provided by DMParameters.
    
    @rtype          :      float
    @return         :      flux.
    """
    x = float(Enu)/float(DMm)
    if x == 0.0 :
        return 0.0
    else : 
        w = np.log10(x)
        pp = p[ch[c]][DM[int(DMm/pc.GeV)]]
        return pp[1]*(1.0+pp[2]*w+pp[3]*w**2+pp[4]*w**3+pp[5]*w**4+pp[6]*w**5)*(1.0-x)**pp[7]+pp[8]*x**pp[9]*(1.0-x)**pp[10]
                        
##END ITALIAN WAY
##BEGIN SWEDISH WAY

def DMSweFlux(Enu,neuflavor,ch,DMm,AtEarth=True):
    """ Gets dN/dz(z) using de swedish data set.
    @type  Enu          :      float
    @param Enu          :      neutrino energy [GeV]
    @type  neuflavor    :      integer
    @param neuflavor    :      flavor of the neutrino (0 : nu_e, 1 : anu_e, 2 : nu_mu 0, 3 : anu_mu, 4 : nu_tau, 5 : anu_tau)
    @type  ch           :      string
    @param ch           :      annihilation channel
    @type  DMm          :      float
    @param DMm          :      dark matter mass [GeV]
    
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
        #z = np.arange(0.0025,0.9975,0.005)
        z=np.linspace(0.0,1.,200)
        h,dat = gt.hreadfilev2(file)
        xneuflavor=neuflavor
        if AtEarth:
            xneuflavor+=12
        dn_dz = dat[xneuflavor]
        #print dat
        
        PC.act_channel, PC.act_DM_mass, PC.act_neuflavor,PC.flag_inter = ch,DMm,neuflavor,False
        
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
        
    def PDF(self,Enu):
        """ Calculates dN/dz where z = E_nu/X_mass using the swedish data.
    
        @type  Enu          :      float
        @param Enu          :      neutrino energy [GeV]
    
        @rtype              :      float
        @return             :      dn/dz [dimensionless]
        """
        return DMSweFlux(Enu,self.flavor,self.ch,self.DMm)#/self.DMm

##END SWEDISH WAY

## DARK MATTER ANNIHILATION RATE ## GOULD ##

def HelmFormFactor(q,A,param):
    """ Calculates Helm's nuclcear form factor. 
    
    Ref : Spin-independent elastic WIMP scattering and the DAMA annual modulation signal. M. Fairbairn.
    arXiv : 0808.0704
    Ec. between ec. (1) and (2).
    
    @type  q            :      float
    @param q            :      transfer momentum [eV]
    @type  A            :      float
    @param A            :      mass number
    @type  param        :      physicsconstants
    @param param        :      set of physical parameters to be used.

    
    @rtype              :      float
    @return             :      Helm's form factor
    """
    s = 1.0*param.fermi
    R = 1.2*A**(1/3)*param.fermi
    r = np.sqrt(R**2-5.0*s**2)
    
    return 3.0*np.exp(-((q*s)**2)/2)*(np.sin(q*r)-q*r*np.cos(q*r))/(q*r)**3
    
def Woods_SaxonFormFactor(q,A,param):
    """ Calculates Wood-Saxong's nuclcear form factor. 
    
    Ref : SUPERSYMMETRIC DARK MATTER, G. JUNGMAN.
    Ec. 7.31 p. 271.
    
    @type  q            :      float
    @param q            :      transfer momentum [eV]
    @type  A            :      float
    @param A            :      mass number
    @type  param        :      physicsconstants
    @param param        :      set of physical parameters to be used.

    
    @rtype              :      float
    @return             :      Woods-Saxon's form factor
    """
    s = 1.0*param.fermi
    R = 1.2*A**(1/3)*param.fermi
    r = np.sqrt(R**2-5.0*s**2)
    
    return (3.0*spe.sph_jn(q*r)/(q*r))*np.exp(-(q*r)**2)
    
def JungmanFormFactor(element_id,DM_mass,param):
    """ Calculates Jungman's form factor according to ec. 9.24 p.301
    
    Ref : SUPERSYMMETRIC DARK MATTER, G. JUNGMAN.
    Ec. 9.24 p. 301.
    
    @type  i            :      int
    @param i            :      element array number.
    @type  DM_mass      :      float
    @param DM_mass      :      DM mass [eV].
    @type  param        :      physicsconstants
    @param param        :      set of physical parameters to be used.

    
    @rtype              :      float
    @return             :      Jungman's form factor
    """
    # USING Table 10, p. 301 from Jungman Rev.
    elements    = {0:'H',1:'He',2:'C',3:'N',4:'O',5:'Ne',6 :'Mg',7:'Si',8:'S',9:'Fe'}   # names and order
    m_e_GeV     = [1.0,18.2,61.6,75.2,75.2,75.2,71.7,71.7,57.0,29.3]                    # [GeV]
    m_e         = [x*param.GeV for x in m_e_GeV]
    F_inf       = [1.0,0.986,0.788,0.613,0.613,0.613,0.281,0.281,0.101,0.00677]
    alpha       = [1.0,1.58,2.69,2.69,2.69,2.69,2.97,2.97,3.1,3.36]
    
    DM_mass_GeV = DM_mass/param.GeV
    
    i = element_id
    if i == 0:                  #Ref : SUPERSYMMETRIC DARK MATTER, G. JUNGMAN. around 9.24 F_H=1
         return 1
    else:
         return F_inf[i]+(1.0-F_inf[i])*np.exp(-(np.log(DM_mass_GeV)/np.log(m_e_GeV[i]))**alpha[i])
    
def GouldAuxF1(A,a,b,eta):
    """ Calculates part of the DM annihilation as estimated by Gould. This is
    the first termn in {}-brackets on eq. 2.18.
    
    Ref : Cosmological density of WIMPs from solar and terrestrial annihilations. A. Gould.
    Astrophysical Journal, 388:338-344, 1992 April.
    
    @type  A      :      float
    @param A      :      Gould's parameter A.
    @type  a      :      float
    @param a      :      Gould's parameter a.
    @type  b      :      float
    @param b      :      Gould's parameter b.
    @type  eta    :      float
    @param eta    :      Gould's parameter eta.
   
    @rtype        :      float
    @return       :      part of Gould's annihilation formulae.
    """
    A_hat   = A*np.sqrt(1.0+a)
    eta_hat = eta/np.sqrt(1.0+a)
    
    A_hat_plus = A_hat + eta_hat
    A_hat_minus= A_hat - eta_hat
    
    return (A_hat_plus*A_hat_minus-0.5-(1+a)/(a-b))*(spe.erf(A_hat_plus)-spe.erf(A_hat_minus))+(1.0/np.sqrt(np.pi))*(A_hat_minus*np.exp(-A_hat_plus**2)-A_hat_plus*np.exp(-A_hat_minus**2))
    
def GouldAuxF2(A,a,b,eta):
    """ Calculates part of the DM annihilation as estimated by Gould. This is
    the first termn in {}-brackets on eq. 2.18.
    
    Ref : Cosmological density of WIMPs from solar and terrestrial annihilations. A. Gould.
    Astrophysical Journal, 388:338-344, 1992 April.
    
    @type  A      :      float
    @param A      :      Gould's parameter A.
    @type  a      :      float
    @param a      :      Gould's parameter a.
    @type  b      :      float
    @param b      :      Gould's parameter b.
    @type  eta    :      float
    @param eta    :      Gould's parameter eta.
   
    @rtype        :      float
    @return       :      part of Gould's annihilation formulae.
    """
    A_hat   = A*np.sqrt(1.0+b)
    eta_hat = eta/np.sqrt(1.0+b)
    
    A_hat_plus = A_hat + eta_hat
    A_hat_minus= A_hat - eta_hat
    
    return (2.0*spe.erf(eta_hat)-spe.erf(A_hat_plus)+spe.erf(A_hat_minus))*np.exp(-(a-b)*A**2)
    

def DMSunAnnihilationRate(DM_mass,DM_cs,param):
    """ Calculates DM annihilation rate in the Sun using Gould's formula.
    
    Ref : Cosmological density of WIMPs from solar and terrestrial annihilations. A. Gould.
    Astrophysical Journal, 388:338-344, 1992 April.
    
    @type  DM_mass      :      float
    @param DM_mass      :      DM mass [eV]
    @type  DM_cs        :      float
    @param DM_cs        :      DM cross section [eV^-2]
    @type  param        :      physicsconstants
    @param param        :      set of physical parameters to be used.

    
    @rtype              :      float
    @return             :      DM annihiliation rate at the Sun. [eV^-1]
    """
    
    # element data
    ''' 
    elements        = {0:'H1', 1:'He4',2:'He3',3:'C12',4:'C13',5:'N14',6:'N15',7:'O16',8:'O17',9:'O18',10:'Ne',11:'Na',12:'Mg',13:'Al',14:'Si',15:'P',16:'S',17:'Cl',18:'Ar',19:'K',20:'Ca',21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',27:'Co',28:'Ni'}          # AGSS09 composition http://www.ice.csic.es/personal/aldos/Solar_Data.html
    #eps 	    = [0.36230,  0.62176,  9.201e-06, 8.386e-06, 1.907e-06, 4.062e-03, 8.275e-08, 6.399e-03, 3.716e-04, 1.640e-05, 1.483e-03, 3.657e-05, 7.035e-04, 6.203e-05, 7.763e-04, 7.121e-06, 3.869e-04, 5.143e-06, 7.789e-05, 4.015e-06, 6.675e-05, 4.309e-08, 3.323e-06, 3.969e-07, 1.939e-05 ,1.417e-05 ,1.344e-03 ,3.732e-06 ,7.946e-05] 
    eps             = [0.75488,  0.23164,  3.313e-05, 2.366e-03, 2.854e-05, 7.009e-04, 1.716e-06, 5.806e-03 ,2.337e-06, 1.396e-05, 1.278e-03, 3.306e-05, 6.360e-04, 5.609e-05, 7.019e-04, 6.438e-06, 3.498e-04, 4.650e-06, 7.042e-05, 3.630e-06, 6.035e-05, 3.895e-08, 3.005e-06, 3.588e-07, 1.753e-05, 1.281e-05, 1.216e-03 ,3.374e-06, 7.184e-05] #AGSS09 composition
    mass_num        = [1.0,4.0,3.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,20.0,23.0,24.0,27.0,28.0,31.0,32.0,35.0,40.0,39.0,40.0,45.0,48.0,51.0,52.0,55.0,56.0,59.0,58.0]
    mass_gr_per_mol = [1.0079,4.0026,3.0160,12.0000,13.0034,14.0031,15.0001,15.9950, 16.9991,17.9992,20.1797,22.9898,24.3050,26.9815,28.0850,30.9738,32.0600,35.4500,39.9481,39.0983,40.0784,44.9559,47.8671,50.9415,51.9962,54.9380,55.8452,58.9332,58.6934]
    '''
    '''
    
    # Bahcall
    elements        = {0:'H',1:'He4',2:'He3',3:'C12',4:'N14',5:'O'}
    mass_num        = [1.0,4.0,3.0,12.0,14.0,16.0]
    eps             = [0.75830,  0.22905,  1.01e-04 , 2.18e-03  ,6.31e-04 , 5.49e-03]
    mass_gr_per_mol = [1.0079,4.0026,3.0160,12.0000,14.0031,15.9950]
    '''
    
    
    #Jungman
    
    elements        = {0:'H',1:'He',2:'C',3:'N',4:'O',5:'Ne',6 :'Mg',7:'Si',8:'S',9:'Fe'}               # names
    mass_num        = [1.0,4.0,12.0,14.0,16.0,20.0,24.0,28.0,32.0,56.0]                                 # A : mass number
    mass_gr_per_mol = [1.0079,4.0026,12.0107,14.0067,15.9994,20.1797,24.3051,28.0855,32.0655,55.8452]   # gr mol^-1
    #eps             = [0.745,0.238,0.0039,0.0009,0.0085,0.0015,0.0007,0.0008,0.0005,0.0016]             # relative aboundances
    eps             = [0.772,0.209,3.87e-3,9.4e-4,8.55e-3,1.51e-3,7.39e-4,8.13e-4,4.65e-4,1.46e-3]       # relative aboundances in the Sun from Table 8 Jungman p.299
    
    n               = len(elements)
    
    # input data
    mass_eV         = [m*param.gr/param.Na for m in mass_gr_per_mol]
    atom_radius     = [(1.2*np.power(A,1.0/3.0)*param.fermi) for A in mass_num]
    
    energy_element  = [3.0/(2.0*mass_eV[i]*atom_radius[i]**2) for i in range(n)]
    
    #DM_rho          = 0.3*param.GeV/param.cm**3 #[Ref : arXiv 1105.6339]
    DM_rho          = 0.39*param.GeV/param.cm**3 # new profile pdg2018
    sun_mass        = 1.9891*1.0e30*param.kg
    
    eta             = 1.0
    # Gould's velocities definitions
    vel_rot         = 220.0*param.km/param.sec          # local DM rel speed.
    vel_char        = vel_rot/eta
    # Jungman velocities definitions
    vel_bar         = 220.0*param.km/param.sec
    vel_star        = np.sqrt(3.0/2.0)*vel_bar/eta
    
    vel_surface     = 795.0*param.km/param.sec
    vel_center      = 1354.0*param.km/param.sec
    
    # special assumptions
    q_i             = [m for m in mass_eV]
    
    # Gould's auxiliaries variables
    beta_plus   = [(4.0*DM_mass*m_i)/(DM_mass+m_i)**2 for m_i in mass_eV]
    beta_minus  = [(4.0*DM_mass*m_i)/(DM_mass-m_i)**2 for m_i in mass_eV]
    #a           = [DM_mass*vel_char**2/(2.0*E) for E in energy_element]
    a           = [DM_mass*vel_bar**2/(2.0*E) for E in energy_element]
    b           = [a[i]*beta_plus[i] for i in range(n)]
    
    #A_center    = [np.sqrt(b_min*(vel_center/vel_char)**2) for b_min in beta_minus]
    #A_surface   = [np.sqrt(b_min*(vel_surface/vel_char)**2) for b_min in beta_minus]
    
    A_center    = [np.sqrt(b_min*(vel_center/vel_bar)**2) for b_min in beta_minus]
    A_surface   = [np.sqrt(b_min*(vel_surface/vel_bar)**2) for b_min in beta_minus]
    
    eta_up_hat  = [eta/np.sqrt(1.0+aa) for aa in a]
    eta_dn_hat  = [eta/np.sqrt(1.0+bb) for bb in b]
    
    # Calculate the annihilation rate as
    # Gamma_i = T_0 * (T_1 + T_2 + T_3)
    # relating Gould's sigma_0 to DM_cs
    reduce_mass = lambda m1,m2 : m1*m2/(m1+m2)
    #sigma_i = [DM_cs*(mass_num[i]*reduce_mass(DM_mass,mass_eV[i])/reduce_mass(DM_mass,p_mass))**2 for i in range(n)]
    
    #sigma_0 = [DM_cs/(beta_plus[0]*mass_eV[0]*DM_mass*q_i[0]**2) for i in range(n)]
    #sigma_0 = [sigma_[i]/(beta_plus[i]*mass_eV[i]*DM_mass*q_i[i]**2) for i in range(n)]
    
    p_mass = 938.272*param.MeV         # proton mass
    q_i = [(mass_num[i]*reduce_mass(DM_mass,mass_eV[i])/reduce_mass(DM_mass,p_mass))/np.sqrt(beta_plus[i]*mass_eV[i]*DM_mass) for i in range(n)]
    sigma_0 = [DM_cs for i in range(n)]
    
    #T_0 = [(sigma_0[i]*DM_rho*sun_mass*vel_char/(2.0*eta))*(q_i[i]**2*eps[i]/a[i]) for i in range(n)]
    T_0 = [(sigma_0[i]*DM_rho*sun_mass*vel_star/(2.0*eta))*(q_i[i]**2*eps[i]/a[i]) for i in range(n)]
    
    T_1 = [(2.0*np.exp(-a[i]*eta_up_hat[i]**2)*spe.erf(eta_up_hat[i])/np.sqrt(1.0+a[i])) for i in range(n)]
    T_2 = [- np.exp(-a[i]*eta_up_hat[i]**2)/((A_center[i]**2-A_surface[i]**2)*np.power(1.0+a[i],3.0/2.0))*(GouldAuxF1(A_center[i],a[i],b[i],eta)-GouldAuxF1(A_surface[i],a[i],b[i],eta)) for i in range(n)]
    T_3 = [np.exp(-b[i]*eta_dn_hat[i]**2)/((a[i]-b[i])*(A_center[i]**2-A_surface[i]**2)*np.sqrt(1.0+b[i]))*(GouldAuxF2(A_center[i],a[i],b[i],eta)-GouldAuxF2(A_surface[i],a[i],b[i],eta)) for i in range(n)]
    
    ann_rate_sun = [T_0[i] * (T_1[i] + T_2[i] + T_3[i]) for i in range(n)]
    
    #return DMSunAnnihilationJungmanSD(DM_mass,DM_cs,param)
    return np.array(ann_rate_sun)/2.0

##

def DMSunAnnihilationRateHooper(DM_mass,DM_cs,param):
    """ Calculates DM annihilation rate in the Sun using Halzen-Cooper formula.
    
    Ref : The Indirect Search for Dark Matter with IceCube. F. Halzen. D. Hooper.
    arXiv : 0910.4513
    
    @type  DM_mass      :      float
    @param DM_mass      :      DM mass [eV]
    @type  DM_cs        :      float
    @param DM_cs        :      DM cross section [eV^-2]
    @type  param        :      physicsconstants
    @param param        :      set of physical parameters to be used.

    
    @rtype              :      float
    @return             :      DM annihiliation rate at the Sun.
    """
    DM_rho          = 0.3*param.GeV/param.cm**3          # local DM density
    #DM_rho          = 0.39*param.GeV/param.cm**3          # local DM density
    vel_rot         = 270.0*param.km/param.sec           # local DM rel speed.
    
    return 3.35*1.0e20*param.sec**-1*(DM_rho/(0.3*param.GeV/param.cm**3))*((270.0*param.km/param.sec)/vel_rot)**3*(100.0*param.GeV/DM_mass)**2*(DM_cs/(1.0e-6*param.picobarn))
    
def DMSunAnnihilationJungman(DM_mass,DM_cs,param):
    """ Calculates DM annihilation rate in the Sun using Jungman/Kamionkowski aproximated formula.
    
    Ref : From Gustav Wikstrom Thesis and paper by Jungman/Kamionkowski et al.
    Eq. 9.20  pp298
    @type  DM_mass      :      float
    @param DM_mass      :      DM mass [eV]
    @type  DM_cs        :      float
    @param DM_cs        :      DM cross section [eV^-2]
    @type  param        :      physicsconstants
    @param param        :      set of physical parameters to be used.

    
    @rtype              :      float
    @return             :      DM annihiliation rate at the Sun [annh/eV^-1].
    """
    #DM_rho          = 0.3*param.GeV/param.cm**3          # local DM density
    DM_rho          = 0.39*param.GeV/param.cm**3          # local DM density
    vel_local       = 270.0*param.km/param.sec           # local DM rel speed.
    
    p_mass          = 938.272*param.MeV         # proton mass
    
    # element data
    elements        = {0:'H',1:'He',2:'C',3:'N',4:'O',5:'Ne',6 :'Mg',7:'Si',8:'S',9:'Fe'}               # names
    n               = len(elements)
    
    mass_num        = [1.0,4.0,12.0,14.0,16.0,20.0,24.0,28.0,32.0,56.0]                                 # A : mass number
    mass_gr_per_mol = [1.0079,4.0026,12.0107,14.0067,15.9994,20.1797,24.3051,28.0855,32.0655,55.8452]   # gr mol^-1
    mass_eV         = [m*param.gr/param.Na for m in mass_gr_per_mol]

    atom_radius     = [(1.2*np.power(A,1.0/3.0)*param.fermi) for A in mass_num]    
    energy_element  = [3.0/(2.0*mass_eV[i]*atom_radius[i]**2) for i in range(n)]
    
    #eps             = [0.745,0.238,0.0039,0.0009,0.0085,0.0015,0.0007,0.0008,0.0005,0.0016]             # relative aboundances in the Sun
    eps             = [0.772,0.209,3.87e-3,9.4e-4,8.55e-3,1.51e-3,7.39e-4,8.13e-4,4.65e-4,1.46e-3]       # relative aboundances in the Sun from Table 8 Jungman p.299
    phi             = [3.16,3.40,3.23,3.23,3.23,3.23,3.23,3.23,3.23,3.23]                                # relative gravitational potential in the Sun from Table 8 Jungman p.299
    # sun
    sun_mass        = 1.9891*1.0e30*param.kg    # Sun mass
    sun_escape      = 1156.0*param.km*param.sec**-1
    
    # auxiliary function
    reduce_mass = lambda m1,m2 : m1*m2/(m1+m2)
    # from Jungman review. ec. 9.21 and 9.22 p.299.
    A_aux = lambda x : (3.0/2.0)*(x/(x-1.0)**2)*(sun_escape/vel_local)**2
    B_aux = 1.5
    kinematical_supression_factor = lambda x : ((A_aux(x)**B_aux)/(1.0+A_aux(x)**B_aux))**(1.0/B_aux)

    nuclear_term = [(((mass_num[i]*reduce_mass(DM_mass,mass_eV[i])/reduce_mass(DM_mass,p_mass))**2)/(mass_eV[i]/param.GeV))*eps[i]*phi[i]*JungmanFormFactor(i,DM_mass,param)*kinematical_supression_factor(mass_eV[i]/DM_mass) for i in range(n)]    
    
    C_c = (4.8e24*param.sec**-1)*(DM_rho/(0.3*param.GeV/param.cm**3))*((270*param.km*param.sec**-1)/vel_local)*(1.0*param.GeV/DM_mass)*(DM_cs/(1.0e-40*param.cm**2))*sum(nuclear_term)
    
    return C_c/2.0

def DMSunAnnihilationJungmanSD(DM_mass,DM_cs,param):
    """ Calculates DM annihilation rate in the Sun using Jungman/Kamionkowski aproximated formula
    in the case of spin dependand interaction.
    
    Ref : From Gustav Wikstrom Thesis and paper by Jungman/Kamionkowski et al.
    Ec. 9.19 on p. 298
    
    @type  DM_mass      :      float
    @param DM_mass      :      DM mass [eV]
    @type  DM_cs        :      float
    @param DM_cs        :      DM cross section [eV^-2]
    @type  param        :      physicsconstants
    @param param        :      set of physical parameters to be used.

    
    @rtype              :      float
    @return             :      DM annihiliation rate at the Sun [annh/eV^-1].
    """
    DM_rho          = 0.3*param.GeV/param.cm**3          # local DM density
    vel_local       = 270.0*param.km/param.sec           # local DM rel speed.
    '''
        
    # element data
    elements        = {0:'H1', 1:'He4',2:'He3',3:'C12',4:'C13',5:'N14',6:'N15',7:'O16',8:'O17',9:'O18',10:'Ne',11:'Na',12:'Mg',13:'Al',14:'Si',15:'P',16:'S',17:'Cl',18:'Ar',19:'K',20:'Ca',21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',27:'Co',28:'Ni'}          # AGSS09 composition http://www.ice.csic.es/personal/aldos/Solar_Data.html
    #eps             = [0.75488,  0.23164,  3.313e-05, 2.366e-03, 2.854e-05, 7.009e-04, 1.716e-06, 5.806e-03 ,2.337e-06, 1.396e-05, 1.278e-03, 3.306e-05, 6.360e-04, 5.609e-05, 7.019e-04, 6.438e-06, 3.498e-04, 4.650e-06, 7.042e-05, 3.630e-06, 6.035e-05, 3.895e-08, 3.005e-06, 3.588e-07, 1.753e-05, 1.281e-05, 1.216e-03 ,3.374e-06, 7.184e-05] #AGSS09 composition
  
  
    mass_num        = [1.0,4.0,3.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,20.0,23.0,24.0,27.0,28.0,31.0,32.0,35.0,40.0,39.0,40.0,45.0,48.0,51.0,52.0,55.0,56.0,59.0,58.0]
    mass_gr_per_mol = [1.0079,4.0026,3.0160,12.0000,13.0034,14.0031,15.0001,15.9950, 16.9991,17.9992,20.1797,22.9898,24.3050,26.9815,28.0850,30.9738,32.0600,35.4500,39.9481,39.0983,40.0784,44.9559,47.8671,50.9415,51.9962,54.9380,55.8452,58.9332,58.6934]
    n               = len(elements)
    
    '''
    
    # element data
    elements        = {0:'H',1:'He',2:'C',3:'N',4:'O',5:'Ne',6 :'Mg',7:'Si',8:'S',9:'Fe'}               # names
    n               = len(elements)
    
    mass_num        = [1.0,4.0,12.0,14.0,16.0,20.0,24.0,28.0,32.0,56.0]                                 # A : mass number
    mass_gr_per_mol = [1.0079,4.0026,12.0107,14.0067,15.9994,20.1797,24.3051,28.0855,32.0655,55.8452]   # gr mol^-1
    
    
    mass_eV         = [m*param.gr/param.Na for m in mass_gr_per_mol]
    
    # sun
    sun_escape      = 1156.0*param.km*param.sec**-1
    
    # from Jungman review. ec. 9.21 and 9.22 p.299.
    A_aux = lambda x : (3.0/2.0)*(x/(x-1.0)**2)*(sun_escape/vel_local)**2
    B_aux = 1.5
    kinematical_supression_factor = lambda x : ((A_aux(x)**B_aux)/(1.0+A_aux(x)**B_aux))**(1.0/B_aux)
    
    C_c = (1.3e25*param.sec**-1)*(DM_rho/(0.3*param.GeV/param.cm**3))*((270*param.km*param.sec**-1)/vel_local)*(1.0*param.GeV/DM_mass)*(DM_cs/(1.0e-40*param.cm**2))*kinematical_supression_factor(DM_mass/mass_eV[0])
    
    return C_c/2.0

## event estimations
    
def DMneuDet(flavor,Enu,ch,DMm,DMsig,body,param,osc):
    """ Calculates Spectra contained event number at detection at a 1 Mton-year Icecube-like detector.
    
    Using the adiabatic aproximation to estimate the probability of survival
    Assuming standard DM parameters for the annhiliation rate we calculate the
    DM-neu spectra at Earth. Then we multiply by the neu-nucleon cross section
    and estimate the number of nuclei in a Mton ice detector. Finally we
    multiply by the exposure and DM sun-annihilation rate.
    
    Calculation :
    
    dN/dE_nu = T * Gamma_ann_rate * Phi_nu(E)_per_ann * cross-section(E) * N_nucleons * P_survival

    @type   flavor  :   integer
    @param  flavor  :   neutrino flavor (0:e,1:mu,2:tau)
    @type   Enu     :   float
    @param  Enu     :   neutrino energy [eV]
    @type   ch      :   string
    @param  ch      :   dark matter annhilation channel
    @type   DMm     :   float
    @param  DMm     :   dark matter mass [eV]
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    @type   osc     :   boolean
    @param  osc     :   if true will consider neu-osc.
    
    @rtype          :   float
    @return         :   dN/dE|detector at E = Enu
    """
    ##B From Arxiv: 0506298 ec. 21 & 24
    #DM_annihilation_rate_Earth  = 1.0e14*(100*param.GeV/DMm)**2/param.sec   #[annhilations/s]
    #DM_annihilation_rate_Sun    = ((1.0*param.AU)/(param.EARTHRADIUS*param.km))**2*DM_annihilation_rate_Earth
    
    DM_annihilation_rate_Sun = float(np.sum(DMSunAnnihilationRate(DMm,DMsig,param))) # [eV^-1]
    
    ##E
    det_exposure    = 1.0*param.year
    det_mass        = 1.0e9*param.kg
    A               = 18.0*param.gr #[gr/mol]
    n_for_n         = 18.0 # [nucleons per nuclei]
    det_nucleons    = n_for_n*(det_mass/A)*param.Na

    evt = 0.0
    if param.neutype == "neutrino":
        if osc :
            for flv in range(3):
                #p = DMParameters(flv)
                if body.name != "vacuum":
                    print param.neutype
                    evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flv*2,ch,DMm/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,param)*param.cm**2*no.AvgNeuProb_RK_STD(flv,flavor,Enu,param)
                    #evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flv*2,ch,DMm/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,0)*param.cm**2*no.AvgNeuProb_RK_STD(flv,flavor,Enu,param)
                    print param.neutype
                else :
                    PMNS = no.mixmatrix(param)
                    evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flv*2,ch,DMm/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,0)*param.cm**2*no.AvgProbNeuVacuum(flv,flavor,PMNS,param)
        else :
            #p = DMParameters(flavor)
            evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flavor*2,ch,DMm/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,0)*param.cm**2
            #evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,0)*param.cm**2#/ycorr
        return evt
    elif param.neutype == "antineutrino":
        if osc :
            for flv in range(3):
                #p = DMParameters(flv)
                if body.name != "vaccum":
                    evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flv*2+1,ch,DMm/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,1)*param.cm**2*no.AvgNeuProb_RK_STD(flv,flavor,Enu,param)
                    #evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,1)*param.cm**2*no.AvgNeuProb_RK_STD(flv,flavor,Enu,param)
                else :
                    PMNS = no.mixmatrix(param)
                    evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flv*2+1,ch,DMm/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,1)*param.cm**2*no.AvgProbNeuVacuum(flv,flavor,PMNS,param)
                    #evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,1)*param.cm**2*no.AvgNeuProb_RK(flv,flavor,Enu,param)#/ycorr
        else :
            #p = DMParameters(flavor)
            evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flavor*2+1,ch,DMm/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,1)*param.cm**2#/ycorr
            #evt = evt + det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,0)*param.cm**-2/ycorr
        return evt
    else :
        print "Wrong neutrino type."
        quit()
        
def DMneuDetMC(flavor,Enu,ch,DMm,DMsig,body,param,inter_flux = None):
    """ Calculates Spectra contained event number at detection at a 1 Mton-year Icecube-like detector.
    
    MC : using montecarlo convoluted data.
    
    Using the adiabatic aproximation to estimate the probability of survival
    Assuming standard DM parameters for the annhiliation rate we calculate the
    DM-neu spectra at Earth. Then we multiply by the neu-nucleon cross section
    and estimate the number of nuclei in a Mton ice detector. Finally we
    multiply by the exposure and DM sun-annihilation rate.
    
    Calculation :
    
    dN/dE_nu = T * Gamma_ann_rate * Phi_nu(E)_per_ann_convoluted_with_probability * cross-section(E) * N_nucleons

    @type   flavor  :   integer
    @param  flavor  :   neutrino flavor (0:e,1:mu,2:tau)
    @type   Enu     :   float
    @param  Enu     :   neutrino energy [eV]
    @type   ch      :   string
    @param  ch      :   dark matter annhilation channel
    @type   DMm     :   float
    @param  DMm     :   dark matter mass [eV]
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    
    @rtype          :   float
    @return         :   dN/dE|detector at E = Enu
    """
    ##B From Arxiv: 0506298 ec. 21 & 24
    #DM_annihilation_rate_Earth  = 1.0e14*(100*param.GeV/DMm)**2/param.sec   #[annhilations/s]
    #DM_annihilation_rate_Sun    = ((1.0*param.AU)/(param.EARTHRADIUS*param.km))**2*DM_annihilation_rate_Earth
    
    DM_annihilation_rate_Sun = float(np.sum(DMSunAnnihilationRate(DMm,DMsig,param)))
    ##E
    det_exposure    = 1.0*param.year
    det_mass        = 1.0e9*param.kg
    A               = 18.0*param.gr #[gr/mol]
    n_for_n         = 18.0 # [nucleons per nuclei]
    det_nucleons    = n_for_n*(det_mass/A)*param.Na

    if inter_flux == None : 
        inter_flux = DMFNeuFluxMCDet(ch,DMm,DMsig,param)

    if param.neutype == "neutrino":
        #return det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*inter_flux(Enu/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,0)*param.cm**2
        return det_nucleons*det_exposure*inter_flux(Enu/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,0)*param.cm**2
    elif param.neutype == "antineutrino":
        #return det_nucleons*det_exposure*(DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*inter_flux(Enu/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,1)*param.cm**2
        return det_nucleons*det_exposure*inter_flux(Enu/param.GeV)*xs.nuDISxsection_NC_NusigmaInt(Enu/param.GeV,1)*param.cm**2
    else :
        print "Wrong neutrino type."
        quit()        

FlavorIndexToState={}
FlavorIndexToState[0]=no.NeuFlavor(1.0,0.0,0.0)
FlavorIndexToState[1]=no.NeuFlavor(0.0,1.0,0.0)
FlavorIndexToState[2]=no.NeuFlavor(0.0,0.0,1.0)
  
def DMFluxneuDet(flavor,Enu,ch,DMm,DMsig,body,param,osc):
    """ Calculates flux spectra number at detection.
    
    Using the adiabatic aproximation to estimate the probability of survival
    Assuming standard DM parameters for the annhiliation rate we calculate the
    DM-neu spectra at Earth.

    @type   flavor  :   integer
    @param  flavor  :   neutrino flavor (0:e,1:mu,2:tau)
    @type   Enu     :   float
    @param  Enu     :   neutrino energy [eV]
    @type   ch      :   string
    @param  ch      :   dark matter annhilation channel
    @type   DMm     :   float
    @param  DMm     :   dark matter mass [eV]
    @type   body    :   body
    @param  body    :   body with the asociated density profile.
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    @type   osc     :   boolean
    @param  osc     :   if true will consider neu-osc.
    
    @rtype          :   float
    @return         :   dPhi/dE_nu|detector at E = Enu
    """                                                                                                        
    ##B From Arxiv: 0506298 ec. 21 & 24
    #DM_annihilation_rate_Earth  = 1.0e14*(100*param.GeV/DMm)**2/param.sec   #[annhilations/s]
    #DM_annihilation_rate_Sun    = ((1.0*param.AU)/(param.EARTHRADIUS*param.km))**2*DM_annihilation_rate_Earth
    DM_annihilation_rate_Sun = float(np.sum(DMSunAnnihilationRate(DMm,DMsig,param)))# [eV]
    ##E
    
    flux = 0.0
    
    if param.neutype == "neutrino":
        if osc :
            PMNS=no.MixMatrix(param)
            for flv in range(3):
                    osc_prob = no.CompleteDecoherenceOscProb(FlavorIndexToState[flv].vec(param),PMNS,param)[flavor]
                #p = DMParameters(flv)
                #if param.name == "STD":
                    flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flv*2,ch,DMm/param.GeV)*osc_prob
                    #flux = flux + (1.0/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)*no.AvgNeuProb_RK_STD(flv,flavor,Enu,param)
                    #flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)*no.AvgNeuProb_RK_STD(flv,flavor,Enu,param)
                #else :
                #    flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flv*2,ch,DMm/param.GeV)*no.AvgNeuProb_RK(flv,flavor,Enu,param)
                    #flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)*no.AvgNeuProb_RK(flv,flavor,Enu,param)
        else :
                #p = DMParameters(flavor)
                flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flavor*2,ch,DMm/param.GeV)
                #flux = flux + (1.0/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)
                #flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)
        return flux
    elif param.neutype == "antineutrino":
        if osc :
            PMNS=no.MixMatrix(param)
            for flv in range(3):
                    osc_prob = no.CompleteDecoherenceOscProb(FlavorIndexToState[flv].vec(param),PMNS,param)[flavor]
                #p = DMParameters(flv)
                #if param.name == "STD":
                    flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flv*2+1,ch,DMm/param.GeV)*osc_prob
                    #flux = flux + (1.0/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)*no.AvgNeuProb_RK_STD(flv,flavor,Enu,param)
                    #flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)*no.AvgNeuProb_RK_STD(flv,flavor,Enu,param)
                #else :
                #    flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flv*2+1,ch,DMm/param.GeV)*no.AvgNeuProb_RK(flv,flavor,Enu,param)
                    #flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)*no.AvgNeuProb_RK(flv,flavor,Enu,param)
        else :
                #p = DMParameters(flavor)
                flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMSweFlux(Enu/param.GeV,flavor*2+1,ch,DMm/param.GeV)
                #flux = flux + (1.0/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)
                #flux = flux + (DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2))*DMFlux(Enu,DMm,ch,p)
        return flux
    else :
        print "Wrong neutrino type."
        quit()
        
        
def DMFNeuFluxMCDet(ch,DMm,DMsig,param):
    """ Calculates muon neutrino flux at the detector using MC data (i.e. considering CC,NC interactions)
    
    Using MC-events we convolute the DM-flux with the propation effects (oscillations,absorptions). We scale
    everything with distance and annhilation rate.

    @type   ch      :   string
    @param  ch      :   dark matter annhilation channel
    @type   DMm     :   float
    @param  DMm     :   dark matter mass [eV]
    @type   DMsig   :   float
    @param  DMsig   :   dark matter scattering cross sectio [eV^-2]
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    
    @rtype          :   spline
    @return         :   return spline interpolating function
    """
    import os
    # FIX SCALING
    ## include years
    DM_annihilation_rate_Sun = DMSunAnnihilationRate(DMm,DMsig,param)       # [eV]
    #DM_annihilation_rate_Sun = 1.6e21/param.sec
    normalization = np.sum((DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2)))  # [eV^3]
    
    ## BEGIN CREATING BINS ##
    # assuming neutrino binnum = 30
    nu_bin_num  = 30
    point_num   = 1000.0
    Emin        = 1.0
    Emax        = 1000.0
    
    E_nu_list   = gt.LogSpaceEnergies(Emin,Emax,binnum = nu_bin_num)
    E_bin_width = [E_nu_list[i+1]-E_nu_list[i] for i in range(len(E_nu_list)-1)]
    E_nu_hpl    = gt.MidPoint(gt.LogSpaceEnergies(Emin,Emax,binnum = nu_bin_num))    
    E_nu_bin    = [0.0]*nu_bin_num # neutrino bins
    E_anu_bin   = [0.0]*nu_bin_num # antineutrino bins
    E_bin_ratio = E_nu_list[1]/E_nu_list[0]
    ## END CREATING BINS ##
    
    for ineu in range(3):
        ## BEGIN READING DATA FROM MC ## 
        
        MCdatapath = "../data/myMC/trials/legion_ineu_"+str(ineu)+"_"+param.name+"/"
        rparam  = PC.PhysicsConstants()
        
        files = []
        for filename in os.listdir(MCdatapath):
            files.append(filename)
                    
        # load all events
        evt = []
        for filename in files :
            file = open(MCdatapath+filename,'r')
            data = []
            gt.hreadfilev4(file,data,rparam)
            if gt.Compareparams(param,rparam):
                print "Using : "+filename
                for e in data :
                    for ee in e:
                        evt.append(ee)
                    
        #del e,ee,data
        
        ## END READING DATA FROM MC ##
        
        # GET DARK MATTER DISTRIBUTION    
        DM_pdf = DM_distribution(ch,DMm/param.GeV,ineu)
           
        for i,e in enumerate(evt):
            if len(e) > 4:
                neutrino = True
                
                family = e[0]
                try:
                    next_family = evt[i+1]
                    if family == next_family and e[1] != 2 :
                        neutrino = False
                except:
                    pass
                
                E_nu_in  = e[2]
                E_nu_out = e[3]
                i = int(np.log(E_nu_out/E_nu_list[0])/np.log(E_bin_ratio))
                j = int(np.log(E_nu_in/E_nu_list[0])/np.log(E_bin_ratio))
                if neutrino:
                    E_nu_bin[i]  = E_nu_bin[i]  + e[5]*(float(DM_pdf.PDF(E_nu_in)/DM_pdf.DMm)*E_bin_width[j]/(np.log(E_nu_list[i])-np.log(E_nu_list[i-1]))) # change to initial neutrino bin width
                    #E_nu_bin[i]  = E_nu_bin[i]  + e[5]*(float(DM_pdf.PDF(E_nu_in)/DM_pdf.DMm))
                else :
                    E_anu_bin[i] = E_anu_bin[i] + e[5]*(float(DM_pdf.PDF(E_nu_in)/DM_pdf.DMm)*E_bin_width[i])
                    #E_anu_bin[i] = E_anu_bin[i] + e[5]*(float(DM_pdf.PDF(E_nu_in)/DM_pdf.DMm))
        
        #int_weight = integrate.quad(lambda E: PDF.PDF(E)/PDF.DMm,Emin,Emax)[0]
        # rescale
    E_nu_bin = [normalization*x/(point_num) for x in E_nu_bin]
    E_anu_bin = [normalization*x/(point_num) for x in E_anu_bin]        
    
    inter_neu = interpolate.InterpolatedUnivariateSpline(E_nu_hpl,E_nu_bin)
    inter_aneu = interpolate.InterpolatedUnivariateSpline(E_nu_hpl,E_anu_bin)
    
    return [inter_neu, inter_aneu]
    
def DMFNeuFluxMCDetv2(ch,DMm,DMsig,param,use_old_data = False,datapath = "../data/myMC/trials/",crosscheck = True,binnum = 200):
    """ Calculates muon neutrino flux number at the detector using MC data (i.e. considering CC,NC interactions)
    
    Using MC-events we convolute the DM-flux with the propation effects (oscillations,absorptions). We scale
    everything with distance and annhilation rate.

    @type   ch      :   string
    @param  ch      :   dark matter annhilation channel
    @type   DMm     :   float
    @param  DMm     :   dark matter mass [eV]
    @type   DMsig   :   float
    @param  DMsig   :   dark matter scattering cross sectio [eV^-2]
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    
    @rtype          :   spline
    @return         :   return spline interpolating function of flux at detector
    """
    import os
    
    existfile_1 = False
    existfile_2 = False
    
    DMm_GeV = DMm/param.GeV
    filename_neutrino = "DMConvolutedFluxes_DMm"+str(DMm_GeV)+"_GeV_DMsig"+str(DMsig/param.cm**2)+"_cm2_channel"+ch+"_"+param.name+"_neutrino.dat"
    filename_antineutrino = "DMConvolutedFluxes_DMm"+str(DMm_GeV)+"_GeV_DMsig"+str(DMsig/param.cm**2)+"_cm2_channel"+ch+"_"+param.name+"_antineutrino.dat"
    
    files = []
    for filename in os.listdir(datapath):
        if filename == filename_neutrino :
            existfile_1 = True
        if filename == filename_antineutrino :
            existfile_2 = True            
            
    if existfile_1 and existfile_2 and use_old_data:
        file_neutrino = open(datapath+filename_neutrino,'r')
        file_antineutrino = open(datapath+filename_antineutrino,'r')

        E_nu_bin    = []
        E_anu_bin   = []

        gt.hreadfilev4(file_neutrino,E_nu_bin,param)
        gt.hreadfilev4(file_antineutrino,E_anu_bin,param)
        
        E_nu_hpl  = map(lambda x : x[0], E_nu_bin[0])
        E_nu_bin  = map(lambda x : x[1], E_nu_bin[0])
        E_anu_bin = map(lambda x : x[1], E_anu_bin[0])
        
    else:
        ## GENERATING DATA FILES ##
        ## BEGIN CREATING BINS ##
        DM_annihilation_rate_Sun = DMSunAnnihilationRate(DMm,DMsig,param)           # [eV]
        normalization = np.sum((DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2)))  # [eV^3]
              
        DMm_GeV = DMm/param.GeV
        
        #print normalization
        #quit()
        # assuming neutrino binnum = 30
        #if param.name != "STD":
        #    nu_bin_num  = 200
        #    point_num   = 1000.0
        #    Emin        = 1.0
        #    Emax        = 10000.0
        #else :
        #    nu_bin_num  = 30
        #    point_num   = 1000.0
        #    Emin        = 1.0
        #    Emax        = 1000.0
            
        nu_bin_num  = binnum
        point_num   = 1000.0
        Emin        = 1.0
        Emax        = 10000.0
        
        E_nu_list   = gt.LogSpaceEnergies(Emin,Emax,binnum = nu_bin_num)
        E_bin_width = [E_nu_list[i+1]-E_nu_list[i] for i in range(len(E_nu_list)-1)]
        E_nu_hpl    = gt.MidPoint(gt.LogSpaceEnergies(Emin,Emax,binnum = nu_bin_num))    
        E_nu_bin    = [0.0]*nu_bin_num # neutrino bins
        E_anu_bin   = [0.0]*nu_bin_num # antineutrino bins
        E_bin_ratio = E_nu_hpl[1]/E_nu_hpl[0]
        
        log_E_bin_ratio = np.log(E_bin_ratio)
        E_nu_list_0 = E_nu_hpl[0]
        

        ## END CREATING BINS ##
        
        # BEGIN GET DARK MATTER DISTRIBUTION
        #flavor of the neutrino (0 : nu_e, 1 : anu_e, 2 : nu_mu 0, 3 : anu_mu, 4 : nu_tau, 5 : anu_tau)
        DM_pdf = []
        for neuneu in range(6):
            DM_pdf.append(DM_distribution(ch,DMm/param.GeV,neuneu))
                
        # calculate DM distribution arrays for futher use
        DM_pdf_table = []
        for neuneu in range(6):
            DM_pdf_table.append(map(lambda EE : float(DM_pdf[neuneu].PDF(EE)/DMm_GeV),E_nu_hpl))
            
        #print [DM_pdf_table[i][120] for i in range(6)]
        #quit()
            
        # END GET DARK MATTER DISTRIBUTION        
    
        for ineu in range(3):
            ## BEGIN READING DATA FROM MC ## 
            # changed for no osc
            #if onlyosc :
            #    MCdatapath = datapath+"legion_ineu_"+str(ineu)+"_"+param.name+"_noosc/"
            #else :
            #    MCdatapath = datapath+"legion_ineu_"+str(ineu)+"_"+param.name+"/"
            MCdatapath = datapath+"legion_ineu_"+str(ineu)+"_"+param.name+"/"
            rparam  = PC.PhysicsConstants()
            
            files = []
            for filename in os.listdir(MCdatapath):
                files.append(filename)
                        
            # load all events for this ineu
            evt = []
            for filename in files :
                print MCdatapath
                print filename
                file = open(MCdatapath+filename,'r')
                data = []
                gt.hreadfilev4(file,data,rparam)
                if crosscheck :
                    if gt.Compareparams(param,rparam):
                        print "Using : "+filename
                        for e in data :
                            for ee in e:
                                evt.append(ee)
                else :
                    print "Using : "+filename
                    for e in data :
                        for ee in e:
                            evt.append(ee)                        
                        
            #del e,ee,data
            
            ## END READING DATA FROM MC ##
            ## BEGIN PROCESS DATA
            E_nu_bin, E_anu_bin = ono.DMEvtProcess(E_nu_bin,E_anu_bin,E_nu_list,E_bin_width,evt,DM_pdf_table,E_nu_list_0,log_E_bin_ratio,DMm_GeV)
            ## THIS CODE WAS DONE FOR OLD DATA COMPATIBILITY ## 
            #if param.name != "STD":
            #    ono.DMEvtProcess(E_nu_bin,E_anu_bin,E_nu_list,E_bin_width,evt,DM_pdf_table,E_nu_list_0,log_E_bin_ratio,DMm_GeV)
            #else :
            #    if param.neutype == "neutrino":
            #        ono.DMEvtProcess_old(E_nu_bin,E_anu_bin,E_nu_list,E_bin_width,evt,DM_pdf_table,E_nu_list_0,log_E_bin_ratio,DMm_GeV,0)
            #    elif param.neutype == "antineutrino":
            #        ono.DMEvtProcess_old(E_nu_bin,E_anu_bin,E_nu_list,E_bin_width,evt,DM_pdf_table,E_nu_list_0,log_E_bin_ratio,DMm_GeV,1)                
            ## END PROCESS DATA
    
        E_nu_bin  = [normalization*x/(point_num) for x in E_nu_bin]
        E_anu_bin = [normalization*x/(point_num) for x in E_anu_bin]
        
        ## begin saving ##
        EE_nu_bin  = []
        EE_anu_bin = []
        
        for i,EE in enumerate(E_nu_hpl):
            EE_nu_bin.append([EE,E_nu_bin[i]])
            EE_anu_bin.append([EE,E_anu_bin[i]])
        
        filename = "DMConvolutedFluxes_DMm"+str(DMm_GeV)+"_GeV_DMsig"+str(DMsig/param.cm**2)+"_cm2_channel"+ch+"_"+param.name+"_neutrino.dat"
        file = open(datapath+filename,'w')
        gt.hwritefile(file,EE_nu_bin,param)
        file.close()
        
        filename = "DMConvolutedFluxes_DMm"+str(DMm_GeV)+"_GeV_DMsig"+str(DMsig/param.cm**2)+"_cm2_channel"+ch+"_"+param.name+"_antineutrino.dat"
        file = open(datapath+filename,'w')
        gt.hwritefile(file,EE_anu_bin,param)
        file.close()
    
    ## end saving ##
    
    inter_neu = interpolate.interp1d(E_nu_hpl,E_nu_bin)
    inter_aneu = interpolate.interp1d(E_nu_hpl,E_anu_bin)
    
    return [inter_neu, inter_aneu]    
    
def DMOscProbabilitiesMC(param,onlyosc = False,datapath = "../data/myMC/trials/",crosscheck = True):
    """ Calculates probabilities from MC

    @type   ch      :   string
    @param  ch      :   dark matter annhilation channel
    @type   DMm     :   float
    @param  DMm     :   dark matter mass [eV]
    @type   DMsig   :   float
    @param  DMsig   :   dark matter scattering cross sectio [eV^-2]
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    
    @rtype          :   spline
    @return         :   return spline interpolating function of flux at detector
    """
    import os
    ## BEGIN CREATING BINS ##
    normalization = 1.0
    # assuming neutrino binnum = 30
    nu_bin_num  = 30
    point_num   = 1000.0
    Emin        = 1.0
    Emax        = 1000.0
    
    E_nu_list   = gt.LogSpaceEnergies(Emin,Emax,binnum = nu_bin_num)
    E_bin_width = [E_nu_list[i+1]-E_nu_list[i] for i in range(len(E_nu_list)-1)]
    E_nu_hpl    = gt.MidPoint(gt.LogSpaceEnergies(Emin,Emax,binnum = nu_bin_num))    
    E_nu_bin    = [0.0]*nu_bin_num # neutrino bins
    E_anu_bin   = [0.0]*nu_bin_num # antineutrino bins
    E_bin_ratio = E_nu_list[1]/E_nu_list[0]
    ## END CREATING BINS ##
    
    for ineu in range(3):
        ## BEGIN READING DATA FROM MC ## 
        # changed for no osc
        
        if onlyosc :
            MCdatapath = datapath+"legion_ineu_"+str(ineu)+"_"+param.name+"_noosc/"
        else :
            MCdatapath = datapath+"legion_ineu_"+str(ineu)+"_"+param.name+"/"
            
        rparam  = PC.PhysicsConstants()
        
        files = []
        for filename in os.listdir(MCdatapath):
            files.append(filename)
                    
        # load all events
        evt = []
        
        for filename in files :
            #print MCdatapath
            print filename
            file = open(MCdatapath+filename,'r')
            data = []
            gt.hreadfilev4(file,data,rparam)
            if crosscheck : 
                if not gt.Compareparams(param,rparam):                
                    print "Using : "+filename
                    for e in data :
                        for ee in e:
                            evt.append(ee)
            else :
                #print "Using : "+filename
                for e in data :
                    for ee in e:
                        evt.append(ee)
                    
        del e,ee,data
        
        ## END READING DATA FROM MC ##
        
        # GET DARK MATTER DISTRIBUTION
        #flavor of the neutrino (0 : nu_e, 1 : anu_e, 2 : nu_mu 0, 3 : anu_mu, 4 : nu_tau, 5 : anu_tau)
        #DM_pdf = []
        #for neutype in range(6):
        #    DM_pdf.append(DM_distribution(ch,DMm/param.GeV,neutype))
           
        for i,e in enumerate(evt):
            #print e
            if len(e) > 4:
                neutrino = True
                
                family = e[0]
                try:
                    next_family = evt[i+1]
                    # cross check this
                    if family == next_family and e[1] != 2 :
                        neutrino = False
                except:
                    pass
                
                iineu = e[1]
                E_nu_in  = e[2]
                E_nu_out = e[3]
                
                j = int(np.log(E_nu_in/E_nu_list[0])/np.log(E_bin_ratio))
                i = int(np.log(E_nu_out/E_nu_list[0])/np.log(E_bin_ratio))
                
                if neutrino:
                    #DMpdf = DM_pdf[int(iineu*2)]
                    #print e
                    if iineu == ineu :
                        E_nu_bin[i]  = E_nu_bin[i]  + e[5]#*(float(DMpdf.PDF(E_nu_in)/DMpdf.DMm))
                    else :
                        pass
                    #E_nu_bin[i]  = E_nu_bin[i]  + e[5]*(float(DMpdf.PDF(E_nu_in)/DMpdf.DMm)*E_bin_width[j])
                else :
                    #DMpdf = DM_pdf[int(iineu*2+1)]
                    #print e
                    E_anu_bin[i] = E_anu_bin[i] + e[5]#*(float(DMpdf.PDF(E_nu_in)/DMpdf.DMm))
                    #E_anu_bin[i] = E_anu_bin[i] + e[5]*(float(DMpdf.PDF(E_nu_in)/DMpdf.DMm)*E_bin_width[j])

    E_nu_bin  = [normalization*x/(1*point_num) for x in E_nu_bin]
    E_anu_bin = [normalization*x/(1*point_num) for x in E_anu_bin]    
    
    inter_neu  = interpolate.InterpolatedUnivariateSpline(E_nu_hpl,E_nu_bin)
    inter_aneu = interpolate.InterpolatedUnivariateSpline(E_nu_hpl,E_anu_bin)
        
    return [inter_neu, inter_aneu]
    
def DMNeuFluxDetNoInt(ch,DMm,DMsig,param,onlyosc = False,datapath = "../data/myMC/trials/"):
    """ Calculates muon neutrino flux number at the detector using analitic oscillation probabilities calculations
    and considering absorbtion effects but neglecting regeneration or NC scattering.

    @type   ch      :   string
    @param  ch      :   dark matter annhilation channel
    @type   DMm     :   float
    @param  DMm     :   dark matter mass [eV]
    @type   DMsig   :   float
    @param  DMsig   :   dark matter scattering cross sectio [eV^-2]
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    @type   onlyosc :   boolean
    @param  onlyosc :   toogle absorption ON/OFF
    @type   datapath:   string
    @param  datapath:   path for the calculated probabilities
    
    @rtype          :   interpolator
    @return         :   return linear interpolating function of pseudo-flux at detector
    """
    DM_annihilation_rate_Sun = DMSunAnnihilationRate(DMm,DMsig,param)           # [eV]
    normalization = np.sum((DM_annihilation_rate_Sun/(4.0*np.pi*param.AU**2)))  # [eV^3]
    
    DM_pdf = []
    for neutype in range(6):
        DM_pdf.append(DM_distribution(ch,DMm/param.GeV,neutype))
    
    E_nu = gt.LogSpaceEnergies(1.0*param.GeV,DMm,binnum = 200)
    
    nu_mu_flux  = []
    anu_mu_flux = []
    
    param.neutype = "neutrino"
    
    for E in E_nu :
        flux = 0.0
        #for ineu in range(param.numneu):
        for ineu in range(3):            
            if onlyosc : 
                PROB = no.AvgNeuProb_RK_STD(ineu,1,E,param,datapath = datapath)
            else :
                print no.AvgNeuProb_RK_STD(ineu,param.numneu + 1,E,param,datapath = datapath)
                print no.NeuSunAbsorptionProbability(E,param)
                PROB = no.NeuSunAbsorptionProbability(E,param)*no.AvgNeuProb_RK_STD(ineu,1,E,param,datapath = datapath)
            DMDIST = (DM_pdf[2*ineu].PDF(E/param.GeV)/DMm)
            XSEC   = xs.nuDISxsection_NC_NusigmaInt(E/param.GeV,0)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)*ice.MuonRange(ice.MuonEnergy(E,0,param),param)*param.meter
            
            flux = flux + normalization*PROB*DMDIST*XSEC
        nu_mu_flux.append(flux)
            
    #print "neutrino success"
    
    param.neutype = "antineutrino"            
            
    for E in E_nu :
        flux = 0.0
        #for ineu in range(param.numneu):
        for ineu in range(3):            
            if onlyosc : 
                PROB = no.AvgNeuProb_RK_STD(ineu,1,E,param,datapath = datapath)
            else :
                print no.AvgNeuProb_RK_STD(ineu,param.numneu + 1,E,param,datapath = datapath)
                print no.NeuSunAbsorptionProbability(E,param)
                PROB = no.NeuSunAbsorptionProbability(E,param)*no.AvgNeuProb_RK_STD(ineu,1,E,param,datapath = datapath)
                quit()
            DMDIST = (DM_pdf[2*ineu+1].PDF(E/param.GeV)/DMm)
            XSEC   = xs.nuDISxsection_NC_NusigmaInt(E/param.GeV,1)*param.cm**2*(0.918*param.gr*param.cm**-3)/(939.27*param.MeV)*ice.MuonRange(ice.MuonEnergy(E,1,param),param)*param.meter
            
            flux = flux + normalization*PROB*DMDIST*XSEC
        anu_mu_flux.append(flux)
            
    total_mu_flux = [float(nu_mu_flux[i]) + float(anu_mu_flux[i]) for i in range(len(nu_mu_flux))]
    
    mu_inter = interpolate.interp1d(E_nu,total_mu_flux)
    
    #print "antineutrino success"
    
    return mu_inter
    
#===========================================================================
# DM FLUX CLASSES
#===========================================================================
        
class DMFluxAtDetector:
    
    def __init__(self,ch,DMm,DMsig,param):
        """ Setup dark matter flux.
   
        @type   ch      :   string
        @param  ch      :   annihilation channel
        @type   DMm     :   float
        @param  DMm     :   dark matter mass [eV]
        @type   param   :   physicsconstants
        @param  param   :   set of physical parameters to be used.
        """
        self.ch = ch
        self.DMm = DMm
        self.DMsig = DMsig        
        self.param = param
        
        self.minenergy = 1.0*param.GeV
        self.maxenergy = DMm
        
    def numu_flux(self,E,neutype):
        """ Returns flux nu_mu event number spectra at detection for a given energy. Uses L{DMFluxneuDet}.
   
        @type   Enu     :   float
        @param  Enu     :   neutrino energy [eV]
        @type   neutype :   integer
        @param  neutype :   neutype  == 0 : neutrino, 1 : antineutrino
       
        @rtype          :   float
        @return         :   dPhi/dE_nu|detector at E = Enu [GeV^-1 m^-2 s^-1]
        """
        Sun = bd.Sun()
        if neutype == 0 :
            self.param.neutype = "neutrino"
        elif neutype == 1:
            self.param.neutype = "antineutrino"
        return DMFluxneuDet(1,E,self.ch,self.DMm,self.DMsig,Sun,self.param,True)*(self.param.GeV*self.param.meter**2*self.param.sec/(self.DMm))
        
class DMMCFluxAtDetector:
    
    def __init__(self,ch,DMm,DMsig,param):
        """ Setup dark matter flux.
   
        @type   ch      :   string
        @param  ch      :   annihilation channel
        @type   DMm     :   float
        @param  DMm     :   dark matter mass [eV]
        @type   param   :   physicsconstants
        @param  param   :   set of physical parameters to be used.
        """
        self.ch = ch
        self.DMm = DMm
        self.DMsig = DMsig
        self.param = param
        
        self.minenergy = 1.0*param.GeV
        self.maxenergy = DMm
        
        self.inter = 0
        
    def numu_flux(self,E,neutype):
        """ Returns flux nu_mu event number spectra at detection for a given energy. Uses L{DMFNeuFluxMCDet}.
   
        @type   Enu     :   float
        @param  Enu     :   neutrino energy [eV]
        @type   neutype :   integer
        @param  neutype :   neutype  == 0 : neutrino, 1 : antineutrino
       
        @rtype          :   float
        @return         :   dPhi/dE_nu|detector at E = Enu [m^-2 s^-1]
        """
        if neutype == 0 :
            self.param.neutype = "neutrino"
        elif neutype == 1:
            self.param.neutype = "antineutrino"
        if self.inter == 0 :
            self.inter = DMFNeuFluxMCDet(self.ch,self.DMm,self.DMsig,self.param)
        return self.inter(E/self.param.GeV)*(self.param.meter**2*self.param.sec)        
 
#===========================================================================
# DM ANNHILATION RATE CROSS CHECK PLOTS
#===========================================================================

def PlotDMCaptureRate(param):
    """ Plot dark matter capture rate at the Sun.
   
    @type   param   :   physicsconstants
    @param  param   :   set of physical parameters to be used.
    """
    plotpath = "../plots/"
    
    DM_mass = np.arange(10.0,1e4,10.0)
    DM_cs   = 1.0e-36*param.cm**2
    #print  DMSunAnnihilationRate(DM_mass[0]*param.GeV,DM_cs,param)
    
    cap_rate_H      = map(lambda dMm : DMSunAnnihilationRate(dMm*param.GeV,DM_cs,param)[0]*pc.sec, DM_mass)
    cap_rate_nuclei = map(lambda dMm : np.sum(DMSunAnnihilationRate(dMm*param.GeV,DM_cs,param))*pc.sec, DM_mass)

    #print DMSunAnnihilationJungmanSD(DM_mass[0]*param.GeV,DM_cs,param)

    #cap_rate_H      = map(lambda dMm : DMSunAnnihilationJungmanSD(dMm*param.GeV,DM_cs,param)*pc.sec, DM_mass)
    #cap_rate_nuclei = map(lambda dMm : np.sum(DMSunAnnihilationJungmanSD(dMm*param.GeV,DM_cs,param))*pc.sec, DM_mass)
   
    plt.plot(DM_mass,cap_rate_H,label = "H", color = "red", linestyle = "solid")
    plt.plot(DM_mass,cap_rate_nuclei,label = "Nuclei", color = "green", linestyle = "solid")
   
    plt.ylabel(r"$C_\odot [\mathrm{s}^{-1}]$")
    plt.xlabel(r"$m_\chi [\mathrm{GeV}]$")
    #plt.semilogy()
    
    plt.ylim(1.0e20,1.0e30)
    plt.xlim(10.0,1000)
    plt.yscale('log')
    #plt.xscale('log')
    plt.grid()
   
    plt.legend(loc = "lower left")
    plt.suptitle("DM capture rate at the Sun")
    
    path = "../plots/"
    filename = "PlotDM_capturerate_sun.png"
    plt.savefig(plotpath+filename)    
    
#===============================================================================
# PATCH
#===============================================================================  
  
def PatchOldDataFiles(newdatapath,olddatapath,param_name = "STD"):
    """ Old version of DM-events were generated using different convention. This
    updates to the new convention.
   
    @type   newdatapath   :   string
    @param  newdatapath   :   datapath of new files
    @type   olddatapath   :   string
    @param  olddatapath   :   datapath of old files
    """
    import os
    
    for neutype in ["neutrino","antineutrino"]:
        
        newdatapath_f = newdatapath + param_name + "/" + neutype +"/"
        olddatapath_f = olddatapath + param_name + "/" + neutype +"/"
        
        for ineu in range(3):
            ## BEGIN READING DATA FROM MC ## 
            MCdatapath = newdatapath_f+"legion_ineu_"+str(ineu)+"_"+param_name+"/"
            MCdatapath_old = olddatapath_f+"legion_ineu_"+str(ineu)+"_"+param_name+"/"
            
            rparam  = PC.PhysicsConstants()
    
            files = []
            for filename in os.listdir(MCdatapath_old):
                files.append(filename)
                        
            # load all events for this ineu
            
            for filename in files :
                old_evt = []
                print MCdatapath_old + "--TO--" + MCdatapath
                print filename
                oldfile = open(MCdatapath_old+filename,'r')
                data = []
                gt.hreadfilev4(oldfile,data,rparam)
                oldfile.close()
                print "Using : "+filename
                for e in data :
                    for ee in e:
                        old_evt.append(ee)
                        
                new_evt = []
                
                for i,e in enumerate(old_evt):
                    keep_type = True
                    noprob = False
                    
                    family = int(e[0])
                    iineu  = int(e[1])
                    E_nu_in = e[2]
                    E_nu_out = e[3]
                    
                    
                    if len(e) > 4 :
                        p1 =e[4]
                        p2 =e[5]
                        p3 =e[6]
                    else :
                        noprob = True
                    
                    try:
                        next_family = int(old_evt[i+1][0])
                        
                        if family == next_family and e[1] != 2 :
                            keep_type = False
                    except:
                        pass
                    
                    if keep_type:
                        if neutype == "neutrino":
                            if noprob :
                                new_evt.append([family,iineu,0,E_nu_in,E_nu_out])
                            else :
                                new_evt.append([family,iineu,0,E_nu_in,E_nu_out,p1,p2,p3])
                        else :
                            if noprob : 
                                new_evt.append([family,iineu,1,E_nu_in,E_nu_out])
                            else :
                                new_evt.append([family,iineu,1,E_nu_in,E_nu_out,p1,p2,p3])
                    else :
                        if neutype == "neutrino":
                            if noprob :
                                new_evt.append([family,iineu,0,E_nu_in,E_nu_out])
                            else :
                                new_evt.append([family,iineu,1,E_nu_in,E_nu_out,p1,p2,p3])
                        else :
                            if noprob : 
                                new_evt.append([family,iineu,1,E_nu_in,E_nu_out])
                            else :
                                new_evt.append([family,iineu,0,E_nu_in,E_nu_out,p1,p2,p3])
                
                newfile = open(MCdatapath+filename,'w')
                gt.hwritefile(newfile,new_evt,rparam)
                newfile.close()
            
  
#===============================================================================
# Testing
#===============================================================================
    
if __name__ == '__main__':
    newdatapath = "../data/myMC_v2/"
    olddatapath = "../data/myMC/"
    
    PatchOldDataFiles(newdatapath,olddatapath,param_name = "STD")
    quit()
    
    pc = PC.PhysicsConstants()
    #p = DMParameters(1)
    #print DMFlux(0.67308e2*pc.GeV,100*pc.GeV,'bb',p)
    
    PlotDMCaptureRate(pc)
    #quit()
    ch = 'tautau'
    neuflavor = 1
    DMm = 100
    
    Enu = 10.0*pc.GeV
    DM_mass = 1000.0*pc.GeV
    DM_csHol   = 1.0e-40*pc.cm**2 # prev. -36
    
    #DMFNeuFluxMCDet(Enu,ch,DM_mass,DM_cs,pc)
    #quit()
    #print DMSunAnnihilationRate(DM_mass,DM_cs,pc)[0]*pc.sec
    print np.sum(DMSunAnnihilationRate(DM_mass,DM_cs,pc)*pc.sec)
    print DMSunAnnihilationRateHooper(DM_mass,DM_cs,pc)*pc.sec
    quit()
    
    #dmflux = DMFluxAtDetector('tautau',DMm*pc.GeV,pc)
    #Enu = np.arange(1.0,DMm,0.1)
    ##Emu = map(lambda E : xs.ymuCC(E,0)*E, Enu)
    #DMSpectra = ice.NevtSpectra(dmflux,0,1.0,pc)
    #import matplotlib.pyplot as plt
    ##plt.plot(Emu,DMSpectra)
    #plt.plot(Enu,DMSpectra)
    #plt.semilogx()
    #plt.xlim(3.0,xs.ymuCC(DMm,0)*DMm)
    #plt.xlim(3.0,DMm)
    #plt.savefig("neu_DM_Icecube.png")
    #print DMm
    ##print xs.ymuCC(DMm,0)*DMm
    ##print xs.dymuCCneu(DMm)*DMm+xs.ymuCC(DMm,0)
    #
    ##print dmflux.flux(1.0*pc.GeV,0)
    
    Sun = bd.Sun()
    
    Enu = 40.0
    DMm = 200.0
    #Enu = np.arange(1.0,DMm,0.1)
    #Enu = np.arange(0.1,100,0.1)
    #xsec = map(lambda x : xs.nuDISxsection_NC_NusigmaInt(x,0)/x,Enu)
    #plt.plot(Enu,xsec)
    #plt.semilogx()
    #plt.savefig("xsection.png")
    
    #DM_events = map(lambda EE : DMneuDet(1,EE*pc.GeV,'tautau',DMm*pc.GeV,Sun,pc,True)*(np.log(10)*Enu/DMm),Enu)
    
    print DMneuDet(1,Enu*pc.GeV,ch,DMm*pc.GeV,Sun,pc,True)*(np.log(10)*Enu/DMm)
