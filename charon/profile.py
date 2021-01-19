"""
Author: Q.R. Liu
"""

import numpy as np
from sympy.solvers import solve
from sympy import Symbol


def NFW(r, rs=24.42, rhos=0.184, gamma=1):
    return rhos / ((r / rs) ** gamma * (1.0 + (r / rs)) ** (3 - gamma))


def Einasto(r, rs=28.44, rhos=0.033, alpha=0.17):
    return rhos * np.exp(-(2.0 / alpha) * ((r / rs) ** alpha - 1.0))


def Burkert(r, rs=12.67, rhos=0.712):
    return rhos / ((1.0 + r / rs) * (1.0 + r ** 2 / rs ** 2))


def Isothermal(r, rs=4.38, rhos=1.387):
    return rhos / (1.0 + r ** 2 / rs ** 2)


# Zhao, H. 1996, MNRAS, 278, 488
# Moore: https://arxiv.org/pdf/astro-ph/0402267.pdf
def Zhao(r, rs=30.28, rhos=0.105, alpha=1.0, beta=3.0, gamma=1.16):
    return rhos / (
        (r / rs) ** gamma * (1.0 + r ** alpha / rs ** alpha) ** ((beta - gamma) / alpha)
    )


def r(theta, l, d):
    # theta: angle between the line of sight and axis connecting Earth and GC in radians
    # l    : distance to the location r from Earth
    # d    : distance from Earth to GC
    # return distance between the location and GC
    return np.sqrt(d ** 2 + l ** 2 - 2 * d * l * np.cos(theta))


def l_max(theta, R, d):
    # R    : size of the galaxy
    x = Symbol("x")
    solution = solve((x * x + d * d - R * R) / (2 * x * d) - np.cos(theta), x)
    return np.float(max(solution))


class J:
    def __init__(self, pro, R, d, process, **kwargs):
        self.pro = pro
        self.R = R
        self.d = d
        self.process = process
        self.params = {}
        if len(kwargs.items()) != 0.0:
            for key, value in kwargs.items():
                self.params[key] = value

    def density(self, r):
        return self.pro(r, **self.params)

    # def Jtheta(pro,theta,R,d,process = 'annihilation',**kwargs):
    def Jtheta(self, theta):
        # return jfactor with unit GeV^2/cm^5 for annihilation and GeV/cm^2 for decay
        l = np.linspace(0.0, l_max(theta, self.R, self.d), 1001)
        width = np.diff(l) * 3.0857e21
        center = (l[1:] + l[:-1]) / 2.0
        r_GC = r(theta, center, self.d)
        densities = self.pro(r_GC, **self.params)
        if self.process == "ann":
            return sum(densities * densities * width)
        elif self.process == "decay":
            return sum(densities * width)

    def J(self, theta_min, theta_max):
        cos_min = np.cos(theta_max)
        cos_max = np.cos(theta_min)
        cos_mid = (cos_min + cos_max) / 2.0
        return 2 * np.pi * self.Jtheta(np.arccos(cos_mid)) * (cos_max - cos_min)

    def JIntegral(self, theta_min, theta_max):
        cos_min = np.cos(theta_max)
        cos_max = np.cos(theta_min)
        cos = np.linspace(cos_min, cos_max, 101)
        cos_mid = (cos[1:] + cos[:-1]) / 2.0
        diff = np.diff(cos)
        return (
            2
            * np.pi
            * sum(
                [
                    self.Jtheta(np.arccos(cos_mid[i])) * (diff[i])
                    for i in range(len(diff))
                ]
            )
        )
