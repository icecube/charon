"""
author : Q.R. Liu
"""

import os, glob
import sys
import time
import socket
from copy import deepcopy

import multiprocessing as multip
from functools import partial

import numpy as np
from scipy.interpolate import interp1d, interp2d
import sympy as sym
from sympy import Symbol
from sympy.solvers import solve

# patch change in nuSQuIDS naming convention
try:
            import nuSQUIDSpy as nsq
except ModuleNotFoundError:
            import nuSQuIDS as nsq

from . import physicsconstants as PC


import astropy.units as u
from astropy.time import Time

import h5py


dirpath = os.path.dirname(os.path.realpath(__file__))
pc = PC.PhysicsConstants()

flavor = {
    0: "nu_e",
    1: "nu_e_bar",
    2: "nu_mu",
    3: "nu_mu_bar",
    4: "nu_tau",
    5: "nu_tau_bar",
}
##############################Multiprocessing##################################
def process(func, args, ncpu):
    def wrapper(func, args, i=0, queue=None):
        r"""function wrapper for multiprocessing

        Parameters
        ----------
        func : function
        args : list
          list of function args where every entry is [tuple, kwargs]
        i : int
          job id
        queue : multiprocessing.Queue

        Returns
        -------
        val : list
          optional list of return values when queue = None
        """
        if isinstance(args[0][0], tuple):
            val = [func(*a[0], **a[1]) for a in args]
        else:
            val = [func(a[0], **a[1]) for a in args]

        # return values if no queue specified
        if queue is None:
            return val

        # otherwise stick answers in the queue
        queue.put([i, val])

    # END wrapper()

    # just return answers for single cpu case
    if ncpu == 1:
        return wrapper(func, args)

    try:
        if os.nice(0) == 0:
            os.nice(1)
    except BaseException:
        warning.warn("Can't set niceness on this computer")
    q = multip.Queue()
    p = [
        multip.Process(target=wrapper, args=(func, a, i, q))
        for i, a in enumerate(np.array_split(args, ncpu))
        if a.size > 0
    ]
    [proc.start() for proc in p]
    r = {}
    for proc in p:
        i, v = q.get()
        r[i] = v
    [proc.join() for proc in p]
    v = []
    for i in range(len(r)):
        v += r[i]
    return v


##############################Multiprocessing##################################

#################################Boost#########################################

def boost(E,DMm,mass_v,theta=0.):
    r"""Boost the energy of products of the mediator decay to lab frame 

    Parameters
    ----------
    E        :  float or narray.    
                energy in rest frame of the mediator in GeV
    DMm      :  float 
                mass of DM  in GeV
    mass_v   :  float 
                mass of mediator in GeV
    theta    :  float
                angle between mediator moving direction and SM products moving direction in radian
    Returns
    -------
    E_lab    :  list 
                energy in lab frame with different angles in CM frame
    dcosCM   :  array           
                widths of angles in CM frame
    """ 
    gamma_v       = DMm/mass_v
    beta_v        = np.sqrt(1.-((1./gamma_v)*(1./gamma_v)))
    momenta       = E
    cosCM         = np.linspace(-1.,1.,101)
    dcosCM        = np.diff(cosCM)
    cosCM         = (cosCM[1:]+cosCM[:-1])/2.
    E_lab         = np.array([E*gamma_v*(1.+beta_v*i) for i in cosCM])
    
    return E_lab, dcosCM

#################################Boost#########################################

##############################Sun/Earth Profile################################
def Model(r, body, path=None):
    r"""Density of matter at distance r to the center

    Parameters
    ----------
    r    :  float
        distance to the center
    body :  string
        "Sun" or "Earth"
    path :  string
        Path to a custom model file. The format should match the default Sun or Earth model file.

    Returns
    -------
    rho  :  float
        density of location r
    """
    if body == "Sun":
        if path == None:
            model = np.genfromtxt(dirpath + "/models/struct_b16_agss09.dat")
        else:
            model = np.genfromtxt(path)
        radius = pc.SUNRADIUS
        radius = model[:, 1] * radius
        rho = model[:, 3]
    elif body == "Earth":
        if path == None:
            model = np.genfromtxt(dirpath + "/models/EARTH_MODEL_PREM.dat")
        else:
            model = np.genfromtxt(path)
        radius = pc.EARTHRADIUS
        radius = model[:, 0] * radius
        rho = model[:, 1]
    f = interp1d(radius, rho, fill_value="extrapolate")
    if r > 0.0 and r <= radius[-1]:
        return float(f(r))
    elif r == 0.0:
        return rho[0]
    else:
        return 0.0


##############################Sun/Earth Profile################################

#################################Geometry######################################


def SunZenith(MJD, lat_det):
    r"""Get the zenith angle of the Solar flux from MJD and detector latitude

    Parameters
    ----------
    MJD       :  float
                Modified Julian Date (MJD)
    lat_det     :  float
                latitude of the detctor in degree
    Returns
    -------
    zenith    :  float
                 zenith angle of the flux in radian
    """
    JD = MJD + 2400000.5
    n = JD - 2451545.0
    L = 280.460 + 0.9856474 * n
    g = 357.528 + 0.9856003 * n
    Lambda = L + 1.915 * np.sin(np.deg2rad(g)) + 0.020 * np.sin(np.deg2rad(2 * g))
    number = Lambda // 360
    Lambda = Lambda - number * 360.0
    epsilon = 23.439 - 0.0000004 * n
    delta = np.arcsin(np.sin(np.deg2rad(epsilon)) * np.sin(np.deg2rad(Lambda)))
    t = Time(MJD, format="mjd", scale="utc", location=(0.0 * u.deg, lat_det * u.deg))
    h = np.deg2rad(
        15
        * (
            (
                float(t.iso[11:13])
                + float(t.iso[14:16]) / 60.0
                + float(t.iso[17:19]) / 3600.0
            )
            - 12.0
        )
    )
    cos = np.sin(np.deg2rad(lat_det)) * np.sin(delta) + np.cos(
        np.deg2rad(lat_det)
    ) * np.cos(delta) * np.cos(h)
    return np.arccos(cos)  # zenith


def Distance(zenith):
    r"""Distances

    Parameters
    ----------
    zenith        :  float
                     zenith angle in radian
    Returns
    -------
    d             :  narray
                     array of [total distance, vacuum distance, atmosphere+earth distance]
    """
    AU = pc.AU / pc.km
    r_earth = pc.EARTHRADIUS
    x = Symbol("x")
    solution = solve(
        (x * x + r_earth * r_earth - AU * AU) / (2 * x * r_earth)
        - sym.cos(sym.pi - zenith),
        x,
    )
    d_tot = [f for f in solution if f > 0][0]
    d_earthatm = nsq.EarthAtm.Track(zenith)
    d_earthatm = d_earthatm.GetFinalX() / pc.km
    d_vacuum = d_tot - d_earthatm
    d = np.array([d_tot, d_vacuum, d_earthatm])
    return d


def xini_Sun(b, r, zenith):
    r"""distances for sun secluded case

    Parameters
    ----------
    b             :  float
                     impact parameter in km
    r             :  float
                     distance to the center in km
    zenith        :  float
                     zenith angle in radian
    Returns
    -------
    d             :  narray
                     array of [total distance, vacuum distance, atmosphere+earth distance]
    """
    R = pc.SUNRADIUS
    l_tot_C = np.float(Distance(zenith)[0])  # Sun center to detector
    d = np.sqrt(l_tot_C ** 2 - b ** 2)
    mid_l = np.sqrt(R ** 2 - b ** 2)
    mid_s = np.sqrt(r ** 2 - b ** 2)
    d_earthatm = nsq.EarthAtm.Track(zenith)
    d_earthatm = d_earthatm.GetFinalX() / pc.km
    if r < R:
        if b == 0:
            return [(-mid_s, mid_s), d - mid_l - d_earthatm]
        else:
            return [(mid_l - mid_s, mid_l + mid_s), d - mid_l - d_earthatm]
    else:
        return [(d + mid_s - d_earthatm, d - mid_s - d_earthatm)]


def xini_Earth(zenith, r):
    r"""distances for Earth secluded case

    Parameters
    ----------
    zenith        :  float
                     zenith angle in radian
    r             :  float
                     distance to the center
    Returns
    -------
    d             :  narray
                     return two initial locations [xini_long_track,xini_short_track] (when r is equal to R_earth*sin(zenith), the two are the same)
    """

    R = pc.EARTHRADIUS
    mid_s = np.sqrt(r * r - (R * np.sin(np.pi - zenith)) ** 2)
    mid_l = R * np.cos(np.pi - zenith)
    return [(mid_l - mid_s, mid_l + mid_s)]


#################################Geometry######################################


#################################Flux Packing#########################################


def Pack(ch, DMm, mass_v, process, folder):
    r"""pack secluded fluxes to shape (len(x),len(densities)). The first density must be 0 which corresponds to vacuum case.
    Parameters
    ----------
    ch              : channel
    DMm             : mass of DM in GeV
    mass_v          : mass of mediator in GeV
    folder          : folder of files
    """

    def key(s):
        return float(
            s[: -s[::-1].find("_") - 1][-s[: -s[::-1].find("_") - 1][::-1].find("_") :]
        )

    flux_list = {}
    for i in range(6):
        files = sorted(
            glob.glob(
                folder + "{}_{:.1f}_{:.1f}_{}_*_--{}.dat".format(ch, DMm, mass_v,process, i)
            ),
            key=lambda s: key(s),
        )
        flux_list[i] = []
        for j in files:
            print (j)
            data = np.genfromtxt(j)
            flux_list[i] += list(data[:, 1])
        flux_list[flavor[i]] = np.transpose(
            np.array(flux_list[i]).reshape(len(files), len(data))
        )
    x_data = data[:, 0]
    flux   = np.zeros(
        (len(data), len(files)),
        dtype=[
            ("x", "float"),
            ("rho", "float"),
            ("nu_e", "float"),
            ("nu_e_bar", "float"),
            ("nu_mu", "float"),
            ("nu_mu_bar", "float"),
            ("nu_tau", "float"),
            ("nu_tau_bar", "float"),
        ],
    )
    for i in range(6):
        flux[flavor[i]] = flux_list[flavor[i]]
    flux["x"] = np.repeat(x_data, len(files)).reshape(len(data), len(files))
    flux["rho"] = np.transpose(np.array([key(s) for s in files] * len(data))).reshape(
        len(data), len(files)
    )
    np.save(folder + "{}_{:.1f}_{:.1f}_{}.npy".format(ch, DMm, mass_v, process), flux)
    return flux


#################################Flux Packing#########################################

###############################Initialize Flux#################################
def LoadFlux(ch, DMm, process="ann"):
    r"""Read Initial Flux and return interpolation functions
    Parameters
    ----------
    ch               :  str
                        Standard Model channel
    DMm              :  float
                        mass of DM in GeV
   
    Returns
    -------
    mass, data       :  array
                        array of masses and corresponding data
    """
    p_mass = {
        "dd": 4.7e-3,
        "uu": 2.2e-3,
        "ss": 95e-3,
        "cc": 1.275,
        "tt": 172.76,
        "bb": 4.18,
        "WW": 80.379,
        "tautau": 1.776,
        "mumu": 105e-3,
        "ee": 0.511e-3,
        "gg": 0.0,
        "ZZ": 91.1876,
        "HH": 125.18,
        "nuenue": 0.0,
        "nutaunutau": 0.0,
        "numunumu": 0.0,
        "gammagamma": 0.0,
    }
    
    if process == "ann":
        factor = 1.0
    elif process == "decay":
        factor = 2.0
   
   
    if DMm / factor < p_mass[ch]:
        sys.exit("DM mass {} GeV is below the threshold of {} channel for {} process".format(DMm, ch, 'a decay' if process == 'decay' else 'an annihilation'))

    sigma_value = 0.56    
    if DMm / factor >= 500.0:
        
        data = h5py.File(dirpath + "/data/SpectraEW.hdf5", "r")
        print("Initial Flux Loading: " + f"Smoothed_SpectraEW_{sigma_value}.hdf5")
        data = h5py.File(data_path + f"Smoothed_SpectraEW_{sigma_value}.hdf5", "r")
        mass = data["m"][:]
    elif DMm / factor < 500.0:
        data = h5py.File(dirpath + "/data/Spectra_PYTHIA.hdf5", "r")
        print("Initial Flux Loading: " + f"Smoothed_Spectra_PYTHIA_{sigma_value}.hdf5")
        data = h5py.File(data_path + f"Smoothed_Spectra_PYTHIA_{sigma_value}.hdf5", "r")
        mass = data["m"][:]
        if (p_mass[ch] >= 3.0):
            mass = np.sort(np.append(p_mass[ch], mass))
        mass = np.unique(mass[mass >= max(3.0, p_mass[ch])])
    return mass, data 



def IniFluxFunction(
    ch, DMm, process="ann", wimp_loc="Sun", path=None, secluded=False, mass_v=1000.0
):
    r"""Read Initial Flux and return interpolation functions
    Parameters
    ----------
    ch               :  str
                        Standard Model channel
    DMm              :  float
                        mass of DM in GeV
    process          :  str
                        annihilation 'ann' or decay 'decay'
    wimp_loc         :  str
                        location of wimp annihilation or decy 'Sun', 'Earth', 'Halo'
    secluded         :  bool
                        standard or secluded
    mass_v           :  float
                        mediator_mass in GeV
    path             :  str
                        path to external file of the flux at production
                        
    Returns
    -------
    f        :  list
                list of interpolation functions of the initial flux for each flavor
    """
   
    if not secluded: 
        f = []
        if not path:
            mass, data = LoadFlux(ch,DMm, process=process)
        
            x = data["x"][:] 
            flux_data = data[wimp_loc][ch]         
            
            for k in range(6):
                f.append(
                    interp2d(
                        x,
                        mass,
                        flux_data[:, k, :],
                        bounds_error=False,
                        kind="linear",
                    )
                )
            data.close()
        
        elif os.path.isfile(path):
            data = np.genfromtxt(path)
            print("Initial Flux Loading: " + path)
            for k in range(6):
                f.append(
                interp1d(data[:, 0], data[:, k + 1], fill_value="extrapolate")
                )
        
        elif os.path.isdir(path):
            for k in range(6):
                if (
                    len(
                        glob.glob(
                            path
                            + "/{}_{}_{}_{}_*-{}.dat".format(
                                ch, DMm, wimp_loc, process, k
                            )
                        )
                    )
                    == 0
                ):
                    sys.exit("Initial flux not found")
                else:
                    for files in glob.glob(
                        path
                        + "/{}_{}_{}_{}_*-{}.dat".format(ch, DMm, wimp_loc, process, k)
                    ):
                        data = np.genfromtxt(files)
                        print("Initial Flux Loading: " + files)
                    f.append(
                        interp1d(data[:, 0], data[:, 1], fill_value="extrapolate")
                    )
        return f
    
    else:       
        if not path:
            loc_list = ['Halo','Earth','Sun']
            mass, data = LoadFlux(ch,mass_v,process="decay")
            xin = data['x'][:]
            x   = xin / 2.
            x   = np.append(x,xin[xin > 0.5])
          
            rhos = np.array([0.0,13.0885,148.9])
            
            flux_data = np.zeros((len(x),3),dtype=[
                ("nu_e", "float"),
                ("nu_mu", "float"),
                ("nu_tau", "float"),
                ("nu_e_bar", "float"),
                ("nu_mu_bar", "float"),
                ("nu_tau_bar", "float"),
            ])
             
            
            
            for k in range(6):
                for j in range(len(rhos)):
                    flux_loc = data[loc_list[j]][ch]                 
                    f_std = interp2d(
                            xin,
                            mass,
                            flux_loc[:, k, :],
                            bounds_error=False,
                            kind="linear",
                            )
                    f_fitted  = 2 * f_std( xin, mass_v / 2.)
                    f_fitted[f_fitted < 0] = 0.0
                    f_fitted[-1] = 0.0                 
                    flux_data[flavor[k]][:,j] = np.append(f_fitted,np.zeros(len(xin[xin > 0.5])))
                    
            data.close()
            
        elif os.path.isfile(path):
            flux_data = np.load(path)
            print("Initial Flux Loading: " + path)
            x        = flux_data["x"][:, 0] / mass_v
            rhos     = flux_data["rho"][0]
       
        elif os.path.isdir(path):
            if os.path.isfile(path + "/{}_{:.1f}_{:.1f}_{}.npy".format(ch, DMm, mass_v, process)):
                data = np.load(path + "/{}_{:.1f}_{:.1f}_{}.npy".format(ch, DMm, mass_v, process))
                print(
                    "Initial Flux Loading: "
                    + path
                    + "/{}_{:.1f}_{:.1f}_{}.npy".format(ch, DMm, mass_v, process)
                )
            else:
                flux_data = Pack(ch, DMm, mass_v, process, path)
                print(
                    "Initial Flux Loading: "
                    + path
                    + "/{}_{:.1f}_{:.1f}_{}.npy".format(ch, DMm, mass_v, process)
                )
           
            x        = flux_data["x"][:, 0] / mass_v
            rhos     = flux_data["rho"][0]
            
        f        = []
        f_vacuum = []

        E_boost, dcosCM   = boost(x * mass_v, DMm, mass_v,theta=0.)


        for k in range(6): 
            boosted_grid = np.zeros((len(rhos),len(x)))
            for rho in range(len(rhos)):
                for l in range(len(E_boost[:,0])):
                    fluxboost = np.zeros(len(x))
                    A  = x * DMm / E_boost[l,:]
                    interp = interp1d(E_boost[l,:] / DMm, 2 * np.pi * A * np.transpose(flux_data[flavor[k]][:, rho]) * abs(dcosCM[l]) / (4 * np.pi),fill_value='extrapolate')                 
                    f_fitted = interp(x)
                    f_fitted[f_fitted < 0] = 0.0
                    boosted_grid[rho,:] += f_fitted
      
              
            f.append(
                interp2d(
                    x,
                    rhos[1:],
                    2 * boosted_grid[1:,:],
                    bounds_error=False,
                )
            )
            f_vacuum.append(
                interp1d(
                    x,
                    2 * boosted_grid[0,:],
                    fill_value="extrapolate",
                )
            )    
        return f, f_vacuum
        
            
def IniFlux(
    Enu,
    ch,
    DMm,
    process="ann",
    wimp_loc="Sun",
    secluded=False,
    r=None,
    mass_v=1000.0,
    pathFlux=None,
    pathSunModel=None,
    pathEarthModel=None,
):
    r"""get initial fluxes at production
    Parameters
    ----------
    Enu              :  np.array
                        neutrino energies in GeV
    ch               :  str
                        Standard Model channel
    DMm              :  float
                        mass of DM in GeV
    process          :  str
                        annihilation 'ann' or decay 'decay'
    wimp_loc         :  str
                        location of wimp annihilation or decy 'Sun', 'Earth', 'Halo'
    secluded         :  bool
                        standard or secluded
    r                :  float or list
                        if secluded, the distance picked
    mass_v           :  float
                        mediator_mass in GeV
    pathFlux         :  str
                        path to the flux file or folder (secluded case)
    pathSunModel     :  str
                        path to the Sun profile
    pathEarthModel   :  str
                        path to the Earth profile

    Returns
    -------
    flux     :  numpy.array
                array of initial flux of all flavors
    """

    nuflavor = {
        "nu_e": 1,
        "nu_e_bar": 2,
        "nu_mu": 3,
        "nu_mu_bar": 4,
        "nu_tau": 5,
        "nu_tau_bar": 6,
    }
    nuID = {
        "nu_e": 12,
        "nu_e_bar": -12,
        "nu_mu": 14,
        "nu_mu_bar": -14,
        "nu_tau": 16,
        "nu_tau_bar": -16,
    }
    
    if process == "ann":
        factor = 1.0
    elif process == "decay":
        factor = 2.0
    indices = np.where(Enu <= DMm / factor) #indices of kinematically allowed Enu values       

    if not secluded:
        flux = np.zeros(
            (len(Enu)),
            dtype=[
                ("nu_e", "float"),
                ("nu_mu", "float"),
                ("nu_tau", "float"),
                ("nu_e_bar", "float"),
                ("nu_mu_bar", "float"),
                ("nu_tau_bar", "float"),
            ],
        )
        f = IniFluxFunction(
            ch,
            DMm,
            process=process,
            wimp_loc=wimp_loc,
            path=pathFlux,
            secluded=secluded,
        )

        if pathFlux == None:
            for k in flux.dtype.names:
                flux[k][indices] = factor * f[nuflavor[k] - 1](factor * Enu[indices] / DMm, DMm / factor)
        else:
            for k in flux.dtype.names:
                flux[k][indices] = factor * f[nuflavor[k] - 1](factor * Enu[indices] / DMm)

    elif secluded:
        f, f_vacuum = IniFluxFunction(
            ch,
            DMm,
            process=process,
            wimp_loc=wimp_loc,
            path=pathFlux,
            secluded=secluded,
            mass_v=mass_v,
        )
        if type(r) in [float, int]:
            r = [r]
        flux = np.zeros(
            (len(r), len(Enu)),
            dtype=[
                ("nu_e", "float"),
                ("nu_mu", "float"),
                ("nu_tau", "float"),
                ("nu_e_bar", "float"),
                ("nu_mu_bar", "float"),
                ("nu_tau_bar", "float"),
            ],
        )
        for k in flux.dtype.names:
            if wimp_loc == "Sun":
                rho = np.array([Model(i, wimp_loc, path=pathSunModel) for i in r])
            elif wimp_loc == "Earth":
                rho = np.array([Model(i, wimp_loc, path=pathEarthModel) for i in r])
            if len(rho[rho > 0.0]) > 1:
                flux_in = np.flip(f[nuflavor[k] - 1](Enu / DMm, rho[rho > 0.0]), axis=0)
            elif len(rho[rho > 0.0]) == 1:
                flux_in = f[nuflavor[k] - 1](Enu / DMm, rho[rho > 0.0])
            else:
                flux_in = np.array([])
            if len(rho[rho == 0.0]) > 0.0:
                flux_out = np.array(
                    [f_vacuum[nuflavor[k] - 1](Enu / DMm)] * len(rho[rho == 0.0])
                )
            else:
                flux_out = np.array([])
            flux[k] = np.append(flux_in, flux_out).reshape(len(rho), len(Enu))
    return flux


###############################Initialize Flux#################################

###############################Averaging oscillation###########################

DTYPE = np.float128
CDTYPE = np.complex256


def angles_to_u(theta12, theta13, theta23, delta):
    r"""Convert angular projection of the mixing matrix elements back into the
    mixing matrix elements.
    Parameters
    ----------
    Returns
    ----------
    unitary numpy ndarray of shape (3, 3)
    """
    theta12 = np.deg2rad(theta12)
    theta13 = np.deg2rad(theta13)
    theta23 = np.deg2rad(theta23)
    s12_2 = np.sin(theta12) ** 2
    c13_4 = np.cos(theta13) ** 4
    s23_2 = np.sin(theta23) ** 2
    dcp = CDTYPE(delta * np.pi / 180.0)

    c12_2 = 1.0 - s12_2
    c13_2 = np.sqrt(c13_4)
    s13_2 = 1.0 - c13_2
    c23_2 = 1.0 - s23_2

    t12 = np.arcsin(np.sqrt(s12_2))
    t13 = np.arccos(np.sqrt(c13_2))
    t23 = np.arcsin(np.sqrt(s23_2))

    c12 = np.cos(t12)
    s12 = np.sin(t12)
    c13 = np.cos(t13)
    s13 = np.sin(t13)
    c23 = np.cos(t23)
    s23 = np.sin(t23)

    p1 = np.array([[1, 0, 0], [0, c23, s23], [0, -s23, c23]], dtype=CDTYPE)
    p2 = np.array(
        [
            [c13, 0, s13 * np.exp(-1j * dcp)],
            [0, 1, 0],
            [-s13 * np.exp(1j * dcp), 0, c13],
        ],
        dtype=CDTYPE,
    )
    p3 = np.array([[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]], dtype=CDTYPE)

    u = np.dot(np.dot(p1, p2), p3)
    return u


def u_to_fr(source_fr, matrix):
    """Compute the observed flavour ratio assuming decoherence.
    Parameters
    ----------
    source_fr : list, length = 3
        Source flavour ratio components
    matrix : numpy ndarray, dimension 3
        Mixing matrix
    Returns
    ----------
    Measured flavour ratio
    ----------
    """
    try:
        composition = np.einsum(
            "ai, bi, a -> b",
            np.abs(matrix) ** 2,
            np.abs(matrix) ** 2,
            source_fr,
        )
    except:
        matrix = np.array(matrix, dtype=np.complex256)
        composition = np.einsum(
            "ai, bi, a -> b",
            np.abs(matrix) ** 2,
            np.abs(matrix) ** 2,
            source_fr,
        )
        pass

    return composition


def average(nuSQ, Eout, flavor):
    r"""average out oscillation

    Parameters
    ----------
    nuSQ     :  nuSQuIDs object
                energy in rest frame of the mediator in GeV
    Eout     :  array
                output energy array
    flavor   :  list
                neutrino flavor nu_e:[0,0]; nu_mu:[1,0]; nu_tau:[2,0]; nu_e_bar:[0,1]; nu_mu_bar:[1,1];
                nu_tau_bar:[2,1]
    Returns
    -------
    flux     :  array
                averaged flux of a specific flavor.
    """

    bins = len(Eout)
    Ediff = np.diff(Eout)
    Ediff = np.append(Ediff[0], Ediff)
    flux = []
    for i, E in enumerate(Eout):
        if i == 0.0:
            E_list = np.linspace(E, E + Ediff[0] / 2.0, 10)
        elif i == bins - 1:
            E_list = np.linspace(E - Ediff[-1] / 2.0, E, 10)
        else:
            E_list = np.linspace(E - Ediff[i] / 2.0, E + Ediff[i] / 2.0, 10)
        flux.append(
            np.mean(
                np.array([nuSQ.EvalFlavor(flavor[0], e, flavor[1]) for e in E_list])
            )
        )
    return np.array(flux)


###############################Averaging oscillation###########################


###########################main propagtion functions###########################


def propagate(
    iniflux,
    Ein,
    Eout,
    location_ini,
    location_end,
    theta_12=33.82,
    theta_23=48.6,
    theta_13=8.60,
    delta_m_12=7.39e-5,
    delta_m_13=2.528e-3,
    delta=0.0,
    interactions=True,
    xsec=None,
    secluded=False,
    r=0.0,
    mass_v=0.0,
    zenith=0.0,
    avg=True,
    pathSunModel=None,
    pathEarthModel=None,
):
    r"""propagate neutrino flux

    Parameters
    ----------
    iniflux         :  array
                       initial flux for propagation
    Ein             :  array
                       energies of the initial flux
    Eout            :  array
                       energies of the output flux
    location_ini    :  str
                       location of neutrino production
    location_end    :  str
                       location that the flux propagates to
    pathSunModel    :  str
                       path to the Sun profile
    pathEarthModel  :  str
                       path to the Earth profile
    theta_12        :  float
                       theta_12 in degree
    theta_23        :  float
                       theta_23 in degree
    theta_13        :  float
                       theta_13 in degree
    delta_m_12      :  float
                       delta_m_12 square in eV^2
    delta_m_13      :  float
                       delta_m_13 square in eV^2
    delta           :  float
                       cp phase in degree
    interactions    :  bool
                       whether to included interactions between nodes
    xsec            :  str
                       Use default cross sections if None. Or can use external table of cross section with names    `n(p)_sigma_CC(NC).dat`,`n(p)_dsde_CC(NC).dat`.
    secluded        :  bool
                       standard or secluded
    r               :  float or list
                       if secluded, the distance picked
    mass_v          :   float
                        mediator_mass in GeV
    zenith          :  float
                       zenith angle in radian
    avg             :  bool
                       whether to average out flux
    Returns
    -------
    flux     :  array
                averaged flux of a specific flavor.
    """

    flavor_list = {
        0: "nu_e",
        1: "nu_e_bar",
        2: "nu_mu",
        3: "nu_mu_bar",
        4: "nu_tau",
        5: "nu_tau_bar",
    }

    if xsec == None:
        try:
            xsec = nsq.CrossSectionLibrary()
            xsec.addTarget(nsq.PDGCode.proton,nsq.NeutrinoDISCrossSectionsFromTables(dirpath + "/xsec/nusigma_proton.h5"))
            xsec.addTarget(nsq.PDGCode.neutron,nsq.NeutrinoDISCrossSectionsFromTables(dirpath + "/xsec/nusigma_neutron.h5"))
        except:
            xsec = nsq.NeutrinoDISCrossSectionsFromTables(dirpath + "/xsec/nusigma_")

    else:
        try:
            xsec = nsq.CrossSectionLibrary()
            xsec.addTarget(nsq.PDGCode.proton,nsq.NeutrinoDISCrossSectionsFromTables(xsec + "_proton.h5"))
            xsec.addTarget(nsq.PDGCode.neutron,nsq.NeutrinoDISCrossSectionsFromTables(xsec + "_neutron.h5"))
        except:
            xsec = nsq.NeutrinoDISCrossSectionsFromTables(xsec)
        else:
            sys.exit("Cross section tables cannot be read.")

    xsec.addTarget(nsq.PDGCode.electron,nsq.GlashowResonanceCrossSection())
    
    
    nuSQ = nsq.nuSQUIDS(Ein * pc.GeV, 3, nsq.NeutrinoType.both, interactions, xsec)
    

    nuSQ.Set_IncludeOscillations(True)
    nuSQ.Set_MixingAngle(0, 1, np.deg2rad(theta_12))
    nuSQ.Set_MixingAngle(0, 2, np.deg2rad(theta_13))
    nuSQ.Set_MixingAngle(1, 2, np.deg2rad(theta_23))
    nuSQ.Set_SquareMassDifference(1, delta_m_12)
    nuSQ.Set_SquareMassDifference(2, delta_m_13)
    nuSQ.Set_abs_error(1.0e-11)
    nuSQ.Set_rel_error(1.0e-11)
    nuSQ.Set_CPPhase(0, 2, np.deg2rad(delta))
    nuSQ.Set_ProgressBar(True)
    if interactions:
        nuSQ.Set_TauRegeneration(True)
        nuSQ.Set_GlashowResonance(True)
    
    flux_output = np.zeros(
        len(Eout),
        dtype=[
            ("Energy", "float"),
            ("nu_e", "float"),
            ("nu_mu", "float"),
            ("nu_tau", "float"),
            ("nu_e_bar", "float"),
            ("nu_mu_bar", "float"),
            ("nu_tau_bar", "float"),
            ("zenith", "float"),
        ],
    )
    # Standard case
    if not secluded:
        # from Sun
        if location_ini == "Sun":
            if pathSunModel == None:
                nuSQ.Set_Body(nsq.Sun(dirpath + "/models/struct_b16_agss09.dat"))
            else:
                nuSQ.Set_Body(nsq.Sun(pathSunModel))
            nuSQ.Set_Track(nsq.Sun.Track(0.0, pc.SUNRADIUS * pc.km))
            nuSQ.Set_initial_state(iniflux, nsq.Basis.flavor)
            nuSQ.EvolveState()

            # flux at sun surface
            if location_end == "SunSurface":
                zenith = None

            # flux at 1AU
            elif location_end == "1AU":
                nuSQ.Set_Body(nsq.Vacuum())
                nuSQ.Set_Track(nsq.Vacuum.Track(pc.AU - pc.SUNRADIUS * pc.km))
                nuSQ.EvolveState()

            # flux at detector
            elif location_end == "detector":
                nuSQ.Set_Body(nsq.Vacuum())
                d_vacuum = float(Distance(zenith)[1])
                nuSQ.Set_Track(nsq.Vacuum.Track(d_vacuum * pc.km))
                nuSQ.EvolveState()
                if pathEarthModel == None:
                    nuSQ.Set_Body(nsq.EarthAtm())
                else:
                    nuSQ.Set_Body(nsq.EarthAtm(pathEarthModel))
                nuSQ.Set_Track(nsq.EarthAtm.Track(zenith))
                nuSQ.EvolveState()

        # from Earth
        elif location_ini == "Earth":

            # flux at detector
            if location_end == "detector":
                if pathEarthModel == None:
                    nuSQ.Set_Body(nsq.Earth())
                else:
                    nuSQ.Set_Body(nsq.Earth(pathEarthModel))
                # nuSQ.Set_Body(nsq.Earth())
                nuSQ.Set_Track(
                    nsq.Earth.Track(
                        pc.EARTHRADIUS * pc.km,
                        2 * pc.EARTHRADIUS * pc.km,
                        2 * pc.EARTHRADIUS * pc.km,
                    )
                )
                nuSQ.Set_initial_state(iniflux, nsq.Basis.flavor)
                nuSQ.EvolveState()

        # from Halo
        elif location_ini == "Halo":

            osc_matrix = angles_to_u(theta_12, theta_13, theta_23, delta)
            composition_nu = np.array(
                [
                    u_to_fr(
                        [iniflux[j, 0, 0], iniflux[j, 0, 1], iniflux[j, 0, 2]],
                        osc_matrix,
                    )
                    for j in range(len(Ein))
                ]
            )
            composition_nubar = np.array(
                [
                    u_to_fr(
                        [iniflux[j, 1, 0], iniflux[j, 1, 1], iniflux[j, 1, 2]],
                        osc_matrix,
                    )
                    for j in range(len(Ein))
                ]
            )

            # flux at Earth surface before entering matter
            if location_end == "Earth" or zenith <= np.pi / 2.0:
                for i in range(3):
                    inter_nu = interp1d(Ein, composition_nu[:, i])
                    inter_nubar = interp1d(Ein, composition_nubar[:, i])
                    flux_output[flavor_list[2 * i]] = inter_nu(Eout)
                    flux_output[flavor_list[2 * i + 1]] = inter_nubar(Eout)

                flux_output["Energy"] = Eout
                flux_output["zenith"] = np.array([zenith] * len(Eout))
                return flux_output

            # flux at detector
            elif location_end == "detector":
                if pathEarthModel == None:
                    nuSQ.Set_Body(nsq.Earth())
                else:
                    nuSQ.Set_Body(nsq.Earth(pathEarthModel))
                # nuSQ.Set_Body(nsq.Earth())
                nuSQ.Set_Track(
                    nsq.Earth.Track(2 * abs(np.cos(zenith)) * pc.EARTHRADIUS * pc.km)
                )
                surface_flux = np.zeros((len(composition_nu), 2, 3))
                for i in range(3):
                    surface_flux[:, 0, i] = composition_nu[:, i]
                    surface_flux[:, 1, i] = composition_nubar[:, i]
                nuSQ.Set_initial_state(surface_flux, nsq.Basis.flavor)
                nuSQ.EvolveState()

    # Secluded case
    elif secluded:

        theta = 0.0
        # from Sun
        if location_ini == "Sun":
            b_impact = np.sin(theta) * r
            if r < pc.SUNRADIUS:
                if theta <= np.pi / 2.0:
                    xini = xini_Sun(b_impact, r, zenith)[0][1]
                else:
                    xini = xini_Sun(b_impact, r, zenith)[0][0]

                if b_impact != 0:
                    l = 2 * np.sqrt(pc.SUNRADIUS ** 2 - b_impact ** 2)
                    if pathSunModel == None:
                        nuSQ.Set_Body(
                            nsq.SunASnu(dirpath + "/models/struct_b16_agss09.dat")
                        )
                    else:
                        nuSQ.Set_Body(nsq.SunASnu(pathSunModel))
                    nuSQ.Set_Track(
                        nsq.SunASnu.Track(l * pc.km, xini * pc.km, b_impact * pc.km)
                    )

                elif b_impact == 0:
                    if pathSunModel == None:
                        nuSQ.Set_Body(
                            nsq.Sun(dirpath + "/models/struct_b16_agss09.dat")
                        )
                    else:
                        nuSQ.Set_Body(nsq.Sun(pathSunModel))
                    nuSQ.Set_Track(nsq.Sun.Track(xini * pc.km, pc.SUNRADIUS * pc.km))

                nuSQ.Set_initial_state(iniflux, nsq.Basis.flavor)
                nuSQ.EvolveState()

            # flux at sun surface
            if location_end == "SunSurface":
                pass

            # flux at 1AU
            elif location_end == "1AU":
                if r < pc.SUNRADIUS:
                    nuSQ.Set_Body(nsq.Vacuum())
                    nuSQ.Set_Track(nsq.Vacuum.Track(pc.AU - pc.SUNRADIUS * pc.km))
                elif r >= pc.SUNRADIUS and r < pc.AU / pc.km - pc.EARTHRADIUS:
                    nuSQ.Set_Body(nsq.Vacuum())
                    nuSQ.Set_Track(nsq.Vacuum.Track(pc.AU - r * pc.km))
                    nuSQ.Set_initial_state(iniflux, nsq.Basis.flavor)
                nuSQ.EvolveState()

            # flux at detector
            elif location_end == "detector":
                nuSQ.Set_Body(nsq.Vacuum())
                if r < pc.SUNRADIUS:
                    nuSQ.Set_Track(
                        nsq.Vacuum.Track(xini_Sun(b_impact, r, zenith)[1] * pc.km)
                    )
                elif r >= pc.SUNRADIUS and r < pc.AU / pc.km - pc.EARTHRADIUS:
                    nuSQ.Set_Track(
                        nsq.Vacuum.Track(xini_Sun(b_impact, r, zenith)[0][1] * pc.km)
                    )
                    nuSQ.Set_initial_state(iniflux, nsq.Basis.flavor)
                nuSQ.EvolveState()
                if pathEarthModel == None:
                    nuSQ.Set_Body(nsq.EarthAtm())
                else:
                    nuSQ.Set_Body(nsq.EarthAtm(pathEarthModel))
                nuSQ.Set_Track(nsq.EarthAtm.Track(zenith))
                nuSQ.EvolveState()

        # from the Earth
        elif location_ini == "Earth" and location_end == "detector":
            xini = xini_Earth(np.pi, r)[0][1]
            if pathEarthModel == None:
                nuSQ.Set_Body(nsq.Earth())
            else:
                nuSQ.Set_Body(nsq.Earth(pathEarthModel))
            nuSQ.Set_Track(
                nsq.Earth.Track(
                    xini * pc.km, 2 * pc.EARTHRADIUS * pc.km, 2 * pc.EARTHRADIUS * pc.km
                )
            )
            nuSQ.Set_initial_state(iniflux, nsq.Basis.flavor)
            nuSQ.EvolveState()

    flux_output["Energy"] = Eout
    flux_output["zenith"] = np.array([zenith] * len(Eout))
    Eout = Eout * pc.GeV
    if avg:
        flux_output["nu_e"] = average(nuSQ, Eout, [0, 0])
        flux_output["nu_mu"] = average(nuSQ, Eout, [1, 0])
        flux_output["nu_tau"] = average(nuSQ, Eout, [2, 0])
        flux_output["nu_e_bar"] = average(nuSQ, Eout, [0, 1])
        flux_output["nu_mu_bar"] = average(nuSQ, Eout, [1, 1])
        flux_output["nu_tau_bar"] = average(nuSQ, Eout, [2, 1])
    else:
        flux_output["nu_e"] = np.array([nuSQ.EvalFlavor(0, e, 0) for e in Eout])
        flux_output["nu_mu"] = np.array([nuSQ.EvalFlavor(1, e, 0) for e in Eout])
        flux_output["nu_tau"] = np.array([nuSQ.EvalFlavor(2, e, 0) for e in Eout])
        flux_output["nu_e_bar"] = np.array([nuSQ.EvalFlavor(0, e, 1) for e in Eout])
        flux_output["nu_mu_bar"] = np.array([nuSQ.EvalFlavor(1, e, 1) for e in Eout])
        flux_output["nu_tau_bar"] = np.array([nuSQ.EvalFlavor(2, e, 1) for e in Eout])

    return flux_output


def wrap_propagate(x, params):
    r"""wrap propagation for multiprocessing
    Parameters
    ----------
    x                :  numpy arrays
                        (initial flux, r)
    params           :  dict
                        parameters

    Returns
    -------
    flux     :  numpy.array
                array of output fluxes at each location in r
    """

    iniflux, r = x

    return propagate(
        iniflux,
        params["Ein"],
        params["Eout"],
        params["location_ini"],
        params["location_end"],
        r=r,
        theta_12=params["theta_12"],
        theta_13=params["theta_13"],
        theta_23=params["theta_23"],
        delta_m_12=params["delta_m_12"],
        delta_m_13=params["delta_m_13"],
        delta=params["delta"],
        interactions=params["interactions"],
        xsec=params["xsec"],
        secluded=True,
        mass_v=params["mass_v"],
        zenith=params["zenith"],
        avg=params["avg"],
        pathSunModel=params["pathSunModel"],
        pathEarthModel=params["pathEarthModel"],
    )


def Probability_point(x, decayLength_v):
    r"""decay probability at one point
    Parameters
    ----------
    x                :  float or array
                        location
    decayLength_v    :  float
                        boosted decay length of the mediator
    """
    prob = 1.0 / decayLength_v * np.exp(-x / decayLength_v)
    return prob


def Probability_length(xini, xend, decayLength_v):
    r"""decay probability of one length
    Parameters
    ----------
    xini             :  float or array
                        initial location
    xend             :  float or array
                        end location
    decayLength_v    :  float
                        boosted decay length of the mediator
    """
    prob = np.exp(-xini / decayLength_v) - np.exp(-xend / decayLength_v)
    return prob


###########################main propagtion functions###########################

##############################Generate fluxes##################################
class NuFlux:
    def __init__(
        self,
        ch,
        DMm,
        nodes,
        Emin=0.0,
        Emax=None,
        bins=None,
        logscale=False,
        process="ann",
        theta_12=33.82,
        theta_23=48.6,
        theta_13=8.60,
        delta_m_12=7.39e-5,
        delta_m_13=2.528e-3,
        delta=0.0,
        interactions=True,
        xsec=None,
        pathFlux=None,
        pathSunModel=None,
        pathEarthModel=None,
    ):
        self.flavor_list = {
            0: "nu_e",
            1: "nu_e_bar",
            2: "nu_mu",
            3: "nu_mu_bar",
            4: "nu_tau",
            5: "nu_tau_bar",
        }
        self.DMm = DMm
        self.ch = ch
        self.nodes = nodes
        self.process = process
        self.interactions = interactions
        self.xsec = xsec
        self.pathFlux = pathFlux
        self.pathSunModel = pathSunModel
        self.pathEarthModel = pathEarthModel
        if bins == None:
            self.bins = nodes
        else:
            self.bins = bins

        if Emax == None:
            Emax = DMm

        if logscale is True:
            self.e_vector = np.logspace(np.log10(Emin), np.log10(Emax), nodes)
            self.e_center = np.logspace(np.log10(Emin), np.log10(Emax), bins)
        else:
            self.e_vector = np.linspace(Emin, Emax, nodes)
            self.e_center = np.linspace(Emin, Emax, bins)

        self.params = {}
        self.params["Ein"] = self.e_vector
        self.params["Eout"] = self.e_center
        self.params["theta_12"] = theta_12
        self.params["theta_13"] = theta_13
        self.params["theta_23"] = theta_23
        self.params["delta_m_12"] = delta_m_12
        self.params["delta_m_13"] = delta_m_13
        self.params["delta"] = delta
        self.params["interactions"] = interactions
        self.params["xsec"] = xsec
        self.params["pathSunModel"] = pathSunModel
        self.params["pathEarthModel"] = pathEarthModel

    def iniE(self):
        r"""Initial energies of the spetra for propagating."""
        return self.e_vector

    def iniFlux(self, wimp_loc):
        r"""initial flux
        wimp_loc   :    string_
                        location of neutrino production
        """
        return IniFlux(
            self.params["Ein"],
            self.ch,
            self.DMm,
            process=self.process,
            wimp_loc=wimp_loc,
            pathFlux=self.pathFlux,
        )

    def iniFlux_secluded(self, wimp_loc, mass_v, r):
        r"""2d initial flux of secluded case
        wimp_loc        :   string
                            object of neutrino production
        r               :   list
                            list of locations picked for interpolation
        mass_v          :   float
                            mediator_mass in GeV
        """
        return IniFlux(
            self.params["Ein"],
            self.ch,
            self.DMm,
            process=self.process,
            wimp_loc=wimp_loc,
            pathFlux=self.pathFlux,
            pathSunModel=self.pathSunModel,
            pathEarthModel=self.pathEarthModel,
            secluded=True,
            r=r,
            mass_v=mass_v,
        )

    def ini_flux(self, wimp_loc):
        r"""initial flux in the form that nuSQuIDs can read"""
        iniflux = self.iniFlux(wimp_loc)
        initial_flux = np.zeros((self.nodes, 2, 3))
        for i in range(3):
            initial_flux[:, 0, i] = iniflux[self.flavor_list[i * 2]]
            initial_flux[:, 1, i] = iniflux[self.flavor_list[i * 2 + 1]]
        return initial_flux

    def ini_flux_secluded(self, wimp_loc, mass_v, r):
        r"""initial flux in the form that nuSQuIDs can read"""
        iniflux = self.iniFlux_secluded(wimp_loc, mass_v, r)
        if type(r) == float:
            r = np.array([r])
        initial_flux = np.zeros((len(r), self.nodes, 2, 3))
        for i in range(3):
            initial_flux[:, :, 0, i] = iniflux[self.flavor_list[i * 2]]
            initial_flux[:, :, 1, i] = iniflux[self.flavor_list[i * 2 + 1]]
        return initial_flux

    def Sun(self, location_end, time=57754.0, zenith=None, latitude=-90.0, avg=True):
        r"""propagate flux from the Sun
        Parameters
        ----------
        location_end    :  str
                           where the flux propagate to. "SunSurface" or "1AU" or "detector"
        zenith          :  float
                           zenith angle of the Sun in radians
        time            :  float
                           MJD if zenith is not specified
        latitude        :  float
                           latitude of the detector in degrees
        avg             :  bool
                           whether to average out oscillation
        """
        initial_flux = self.ini_flux("Sun")

        if zenith == None:
            zenith = SunZenith(time, latitude)
        else:
            zenith = zenith
        return propagate(
            iniflux=initial_flux,
            location_ini="Sun",
            location_end=location_end,
            secluded=False,
            zenith=zenith,
            avg=avg,
            **self.params
        )

    def Earth(self, location_end, avg=True):
        r"""propagate flux from the Earth
        Parameters
        ----------
        location_end    :  str
                           where the flux propagate to. "detector"
        avg             :  bool
                           whether to average out oscillation
        """
        initial_flux = self.ini_flux("Earth")
        return propagate(
            iniflux=initial_flux,
            location_ini="Earth",
            location_end=location_end,
            secluded=False,
            avg=avg,
            **self.params
        )

    def Halo(self, location_end, zenith=0.0, avg=True):
        r"""averaged propagated flux from the Halo
        Parameters
        ----------
        location_end    :  str
                           where the flux propagate to. "Earth" or "detector"
        zenith          :  float
                           zenith angle of the flux in radians
        avg             :  bool
                           whether to average out oscillation in Earth
        """
        initial_flux = self.ini_flux("Halo")
        return propagate(
            iniflux=initial_flux,
            location_ini="Halo",
            location_end=location_end,
            secluded=False,
            avg=avg,
            zenith=zenith,
            **self.params
        )

    def SunSecluded(
        self,
        segment,
        lambda_v,
        mass_v,
        location_end,
        ncpu=1,
        time=57754.0,
        zenith=None,
        latitude=-90.0,
        avg=False,
        p=0.9,
        approx_1D=True,
    ):
        r"""propagate flux from the Sun
        Parameters
        ----------
        segment         :  int
                           number of points picked along the line of sight
        lambda_v        :  float
                           decay length in km
        mass_v          :  float
                           mediator_mass in GeV
        location_end    :  str
                           where the flux propagate to. "SunSurface" or "1AU" or "detector"
        ncpu            :  int
                           number of cpus to do multiprocessing
        p               :  float
                           fraction of decays along the line of sight which decides length to compute the flux
        zenith          :  float
                           zenith angle of the Sun in radians
        time            :  float
                           MJD if zenith is not specified
        latitude        :  float
                           latitude of the detector in degrees
        avg             :  bool
                           whether to average out oscillation
        approx_1D       :  bool
                           whether to use 1D approximation
        """
        xmax = -lambda_v * np.log(1 - p)
        r = np.linspace(0.0, xmax, segment + 1)
        r_center = (r[1:] + r[:-1]) / 2.0

        initial_flux = self.ini_flux_secluded("Sun", mass_v, r_center)
        if zenith == None:
            zenith = SunZenith(time, latitude)
        else:
            zenith = zenith

        ss_params = deepcopy(self.params)
        ss_params["location_ini"] = "Sun"
        ss_params["location_end"] = location_end
        ss_params["mass_v"] = mass_v
        ss_params["zenith"] = zenith
        ss_params["avg"] = avg

        pool = multip.Pool(processes=ncpu)
        fluxes = pool.map(
            partial(wrap_propagate, params=ss_params),
            [(initial_flux[i], r_center[i]) for i in range(segment)],
        )
        results = np.zeros(
            (self.bins),
            dtype=[
                ("Energy", "float"),
                ("nu_e", "float"),
                ("nu_mu", "float"),
                ("nu_tau", "float"),
                ("nu_e_bar", "float"),
                ("nu_mu_bar", "float"),
                ("nu_tau_bar", "float"),
                ("zenith", "float"),
            ],
        )
        if approx_1D:
            prob = Probability_length(r[:-1], r[1:], lambda_v)

            for j in results.dtype.names:
                if j != "Energy" or j != "zenith":
                    results[j] = sum(
                        [fluxes[k][j] * prob[k] for k in range(len(r_center))]
                    )
            results["Energy"] = fluxes[0]["Energy"]
            results["zenith"] = fluxes[0]["zenith"]

        return results

    def EarthSecluded(
        self,
        segment,
        lambda_v,
        mass_v,
        location_end,
        zenith=None,
        p=0.9,
        avg=False,
        approx_1D=True,
        ncpu=1,
    ):
        xmax = -lambda_v * np.log(1 - p)
        if approx_1D:
            xmax = min(pc.EARTHRADIUS, xmax)
        r = np.linspace(0.0, xmax, segment + 1)
        r_center = (r[1:] + r[:-1]) / 2.0
        initial_flux = self.ini_flux_secluded("Earth", mass_v, r_center)

        es_params = deepcopy(self.params)
        es_params["location_ini"] = "Earth"
        es_params["location_end"] = location_end
        es_params["mass_v"] = mass_v
        es_params["zenith"] = zenith
        es_params["avg"] = avg

        pool = multip.Pool(processes=ncpu)
        fluxes = pool.map(
            partial(wrap_propagate, params=es_params),
            [(initial_flux[i], r_center[i]) for i in range(segment)],
        )
        results = np.zeros(
            (self.bins),
            dtype=[
                ("Energy", "float"),
                ("nu_e", "float"),
                ("nu_mu", "float"),
                ("nu_tau", "float"),
                ("nu_e_bar", "float"),
                ("nu_mu_bar", "float"),
                ("nu_tau_bar", "float"),
                ("zenith", "float"),
            ],
        )
        if approx_1D:
            prob = Probability_length(r[:-1], r[1:], lambda_v)
            for j in results.dtype.names:
                if j != "Energy" or j != "zenith":
                    results[j] = sum(
                        [fluxes[k][j] * prob[k] for k in range(len(r_center))]
                    )
            results["Energy"] = fluxes[0]["Energy"]
            results["zenith"] = fluxes[0]["zenith"]

        return results


##############################Generate fluxes##################################
