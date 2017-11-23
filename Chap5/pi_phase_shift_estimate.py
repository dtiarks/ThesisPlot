#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 11:25:42 2017

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


#R_b=18e-6

def fitFunc(t, A, phi,C):
    return A*np.sin(2*np.pi*20*t+phi)+C

def susceptibility(Delta_s, Delta_c, gamma_21, Omega_c, ladder=True):
    delta = (Delta_s + (-1 + 2*int(ladder)) * Delta_c)
    return 1j*(gamma_21 - 2j * delta)/(np.abs(Omega_c)**2 + (1 - 2j * Delta_s)*(gamma_21 - 2j * delta))

def susceptibility_off(Delta_s, Delta_c, gamma_21, Omega_c, ladder=True):
    return -1/(2*Delta_s+1j)

def vdW_pot(r, r0):
    return -C6 * (r-r0)**-6

def cond_phase(d_c,om_c,d_o):
    im_chi_vdw = lambda x,y: omega_s/(c)*np.imag(susceptibility(y, d_c - vdW_pot(x, 1e-10)/(hbar*Gamma_e), gamma_21, om_c))
    od_vdw = lambda x: integrate.quad(im_chi_vdw, -L/2, L/2,args=(x,))[0]

    intersection = newton(lambda x: L*omega_s/(c)*np.imag(susceptibility(x, d_c, gamma_21, om_c) - susceptibility(x, d_c, gamma_21, 0))-d_o, -d_c)
#    intersection = newton(lambda x: L*omega_s/(c)*np.imag(susceptibility(x, d_c, gamma_21, om_c)) - od_vdw(x)-d_o, -d_c)

    chi_nb = susceptibility(intersection, d_c, gamma_21, om_c)
    phi_0=omega_s/(2*c) * L * chi_0 * np.real(chi_nb)

    r_chi_vdw = lambda x: np.real(susceptibility(intersection, d_c - vdW_pot(x, 1e-10)/(hbar*Gamma_e), gamma_21, om_c))
    phi_1=omega_s/(2*c) *  chi_0 *integrate.quad(r_chi_vdw, -L/2, L/2)[0]

    d_phi=phi_1-phi_0
    return intersection,d_phi

def cond_trans(d_c,om_c,d_o):
    intersection = newton(lambda x: L*omega_s/(c)*np.imag(susceptibility(x, d_c, gamma_21, om_c) - susceptibility(x, d_c, gamma_21, 0))-d_o, -d_c)

    im_chi_vdw = lambda x: np.imag(susceptibility(intersection, d_c - vdW_pot(x, 1e-10)/(hbar*Gamma_e), gamma_21, om_c))
    t_1=omega_s/c *  chi_0 *integrate.quad(im_chi_vdw, -L/2, L/2)[0]


    return intersection,np.exp(-t_1)


# the parameters, all in units of Gamma_3
Delta_c = 2*np.pi*9.2*10**6/Gamma_e
Delta_s = -2*np.pi*10.*10**6/Gamma_e
ds_off = 2*np.pi*0.0*10**6/Gamma_e

Omega_c = 2*np.pi*10.4*10**6/Gamma_e

#Ga=Omega_c**2/(4*np.abs(Delta_s))
Ga=np.true_divide(Omega_c**2*np.abs(Delta_s),1+np.sqrt((4*np.abs(Delta_s)**2+1)**2-4*np.abs(Delta_s)**2*1))
R_b=1.0*(C6/(hbar*Gamma_e)/Ga)**(1./6.)

print("Blockade Radius: %.2f um"%(R_b*1e6))

od_new = lambda x: omega_s/c*chi_0*L*np.imag(susceptibility(x-ds_off, Delta_c , gamma_21, Omega_c))
ph_new = lambda x: 0.5*omega_s/c*chi_0*L*np.real(susceptibility(x-ds_off, Delta_c , gamma_21, Omega_c))-0.5*omega_s/c*chi_0*L/15.*np.real(susceptibility(x, Delta_c , gamma_21, Omega_c))
#od0_new = lambda x: 0.7*omega_s/(c)* L*chi_0*np.imag(susceptibility(x-ds_off, Delta_c , gamma_21, 0*Omega_c))
od0_new = lambda x: omega_s/c*chi_0*L*np.imag(susceptibility_off(x-0*ds_off, Delta_c , gamma_21, 0*Omega_c))
ph0_new = lambda x: 0.5*omega_s/c*chi_0*L*np.real(susceptibility_off(x-0*ds_off, Delta_c , gamma_21, 0*Omega_c))-0.5*omega_s/c*chi_0*L/15.*np.real(susceptibility_off(x, Delta_c , gamma_21, 0*Omega_c))


od_max=omega_s/c*chi_0*L

print("gamma_21: %.2f"%gamma_21)
print("est. od_max: %.2f"%od_max)

od_max=26.3


chi_nb = susceptibility(Delta_s, Delta_c, gamma_21, Omega_c)
phi_0=0.5*od_max*np.real(chi_nb)

r_chi_vdw = lambda x: np.real(susceptibility(Delta_s, Delta_c - vdW_pot(x, 1e-10)/(hbar*Gamma_e), gamma_21, Omega_c))
phi_1=0.5*od_max*integrate.quad(r_chi_vdw, -L/2, L/2)[0]/L

D_phi_blockaded=phi_1-phi_0

#D_phi=ph_new(Delta_s)-ph0_new(Delta_s)

#D_phi_blockaded = D_phi*2*R_b/L

print("Expected cond. phase shift: %.3f rad"%D_phi_blockaded)

#R_b=15e-6
print("Simple cond. phase shift: %.3f rad"%(6.6*2*R_b/L))


#fig=plt.figure(1, figsize=(8, 9), dpi=80,)
#fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)
#
#ax0 = fig.add_subplot(211)
#ax0.clear()
#ds=np.linspace(4,16,num=n)
#dsf=np.linspace(4,16,num=len(freqFitO))
##dsf2=np.linspace(4,16,num=len(freqFitO))
#ax0.plot(ds,2*np.log(aRef[:,0]/a[:,0]),'ro',label="$\Omega_c=$???")
#ax0.plot(ds,2*np.log(aRef[:,1]/a[:,1]),'bo',label="$\Omega_c=$0")
##ax0.plot(dsf,OfitOff,"b")
##ax0.plot(dsf,OfitOn,"r")
#
#ax0.plot(freqFitP,OfitOff,"b")
#ax0.plot(freqFitP,OfitOn,"r")
#ax0.plot(dsf,od_new(dsf*1e6*2*np.pi/Gamma_e),"g")
#ax0.plot(dsf,od0_new(dsf*1e6*2*np.pi/Gamma_e),"k")
#
#ax0.tick_params(axis="x",which="both",labelbottom="off")
#ax0.set_ylabel("OD")
#trans = ax0.get_xaxis_transform() # x in data untis, y in axes fraction
#ax0.annotate('(a)', xy=(4.3,0.9 ), xycoords=trans)
#
#
#ax1 = fig.add_subplot(212)
#ax1.clear()
#
#
#dsf=np.linspace(4,16,num=len(freqFitP))
#ax1.plot(ds,(p[:,0]-pRef[:,0]),'ro',label="control on")
#ax1.plot(ds,(p[:,1]-pRef[:,1]),'bo',label="control off")
##ax1.plot(dsf,pfitOff,"b")
##ax1.plot(dsf,pfitOn,"r")
#
#ax1.plot(freqFitP,pfitOff,"b")
#ax1.plot(freqFitP,pfitOn,"r")
#ax1.plot(dsf,ph_new(dsf*1e6*2*np.pi/Gamma_e)+0*phase_offset,"g")
#ax1.plot(dsf,ph0_new(dsf*1e6*2*np.pi/Gamma_e)+0*phase_offset,"k")
#
#plt.show()

