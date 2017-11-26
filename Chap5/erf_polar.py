#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 16:05:59 2017

@author: daniel
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import io
from scipy.optimize import newton
from scipy import integrate
from scipy import special
from scipy.optimize import curve_fit
import re


c = 299792458 # m/s, speed of light CODATA 2014
a0 = 0.52917721067e-10 # m, Bohr radius
C6 = 2.3e23 * 4.36e-18 * a0**6 # Jm^6, Van-der-Waals coefficient for the 67s - 69s
hbar = 6.626070040e-34/(2 * np.pi) # Js, Planck constant, CODATA 2014
rho_peak = 1.8e12/1e-6 # peak density in cm^-3/centi^-3
d = 2.534e-29 # Cm, dipole matrix element (D. A. Steck)
Gamma_e = 2*np.pi * 6.065e6 # decay rate (D. A. Steck)
epsilon_0 = 8.854187817e-12 # dielectric constant, CODATA 2014
L = 61e-6 # medium length in m
omega_s = 2*np.pi * 384.23e12 # rad/s, transition frequency
gamma_21 = 2.2/(2*np.pi)
chi_0 = 2*rho_peak*d**2 / (epsilon_0*hbar*Gamma_e) # prefactor of the susceptibility for the cycling transition (|R> polarization)

def pol_komp(t,r_f,r_b, wd,d):
    Ip=1-special.erf((txs-(wd+d))/r_f)*special.erf((txs-d)/r_b)
    return Ip

r_f_plus=0.2
r_b_plus=0.1
w_plus=3.
d_plus=0

r_f_minus=0.2
r_b_minus=0.1
w_minus=3.
d_minus=1.

txs=np.linspace(-0.3,4.5,1000)

I_p=0.8*pol_komp(txs,r_f_plus,r_b_plus, w_plus,d_plus)
I_m=pol_komp(txs,r_f_minus,r_b_minus, w_minus,d_minus)

S3=np.true_divide(I_p-I_m,I_p+I_m+1e-12)

theta=np.arctan(np.clip(1./(S3),-1e12,1e12))


plt.plot(txs, theta)
plt.plot(txs, pol_komp(txs,r_f_plus,r_b_plus, w_plus,d_plus))
plt.plot(txs, pol_komp(txs,r_f_minus,r_b_minus, w_plus,d_minus))

plt.show()
