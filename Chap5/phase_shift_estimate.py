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

phaseLst=[]
phaseLstRef=[]
ampLst=[]
ampLstRef=[]

n=21

c = 299792458 # m/s, speed of light CODATA 2014
a0 = 0.52917721067e-10 # m, Bohr radius
C6 = 2.3e23 * 4.36e-18 * a0**6 # Jm^6, Van-der-Waals coefficient for the 67s - 69s
hbar = 6.626070040e-34/(2 * np.pi) # Js, Planck constant, CODATA 2014
rho_peak = 3.1e12/1e-6 # peak density in cm^-3/centi^-3
d = 2.534e-29 # Cm, dipole matrix element (D. A. Steck)
Gamma_e = 2*np.pi * 5.746e6 # decay rate (D. A. Steck)
epsilon_0 = 8.854187817e-12 # dielectric constant, CODATA 2014
L = 2*24e-6 # medium length in m
omega_s = 2*np.pi * 377e12 # rad/s, transition frequency
gamma_21 = 1.3/5.75
chi_0 = 0.5*2*rho_peak*d**2 / (epsilon_0*hbar*Gamma_e) # prefactor of the susceptibility for the cycling transition (|R> polarization)
chi_0 = 2*c/(2*omega_s*L*9.)

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
#    return intersection,t_1
    
    
# the parameters, all in units of Gamma_3
Delta_c = -2*np.pi*8.*10**6/Gamma_e
Delta_s = 2*np.pi*8.0*10**6/Gamma_e
ds_off = 2*np.pi*0.57*10**6/Gamma_e

Omega_c = 2*np.pi*5.6*10**6/Gamma_e

#Ga=Omega_c**2/(4*Delta_s)
Ga=np.true_divide(Omega_c**2*Delta_s,1+np.sqrt((4*Delta_s**2+1)**2-4*Delta_s**2*1))
R_b=1.0*(C6/(hbar*Gamma_e)/Ga)**(1./6.) 

od_new = lambda x: 7.0*np.imag(susceptibility(x-ds_off, Delta_c , gamma_21, Omega_c))
ph_new = lambda x: 0.5*9.*np.real(susceptibility(x-ds_off, Delta_c , gamma_21, Omega_c))-0.5*7.7*np.real(susceptibility(x+2*np.pi*20*10**6/Gamma_e, Delta_c , gamma_21, Omega_c))
#od0_new = lambda x: 0.7*omega_s/(c)* L*chi_0*np.imag(susceptibility(x-ds_off, Delta_c , gamma_21, 0*Omega_c))
od0_new = lambda x: 8.3*np.imag(susceptibility_off(x-0*ds_off, Delta_c , gamma_21, 0*Omega_c))
ph0_new = lambda x: 0.5*9.*np.real(susceptibility_off(x-0*ds_off, Delta_c , gamma_21, 0*Omega_c))-0.5*7.9*np.real(susceptibility_off(x+2*np.pi*20*10**6/Gamma_e, Delta_c , gamma_21, 0*Omega_c))

for i in range(1,n+1):
   name="Phase_vs_fProbe.tsv_%03d"%i

   histo=np.loadtxt(os.path.join("Giovanni","histo_%s"%name),skiprows=2000)
   histoRef=np.loadtxt(os.path.join("Giovanni","histoRef_%s"%name),skiprows=2000)
   histoInterval=(2530-2000,2580-2000)

   hT=histo[histoInterval[0]:histoInterval[1],0]
   hCrop1 = histo[histoInterval[0]:histoInterval[1],1]
   hCrop2 = histo[histoInterval[0]:histoInterval[1],2]
   
   hTRef=histoRef[histoInterval[0]:histoInterval[1],0]
   hCrop1Ref = histoRef[histoInterval[0]:histoInterval[1],1]
   hCrop2Ref = histoRef[histoInterval[0]:histoInterval[1],2]

   popt0=[0.01, 4,0.002]
   popt1=[0.01, 4,0.002]

   popt0, pcov0=curve_fit(fitFunc, hT, hCrop1, popt0)
   popt1, pcov0=curve_fit(fitFunc, hT, hCrop2, popt1)
   
   popt0Ref, pcov0Ref=curve_fit(fitFunc, hTRef, hCrop1Ref, popt0)
   popt1Ref, pcov0Ref=curve_fit(fitFunc, hTRef, hCrop2Ref, popt1)

   phaseLst.append([popt0[1],popt1[1]])
   ampLst.append([popt0[0],popt1[0]])
   phaseLstRef.append([popt0Ref[1],popt1Ref[1]])
   ampLstRef.append([popt0Ref[0],popt1Ref[0]])
   print("sample "+str(i)+" analyzed")

p=np.array(phaseLst)
a=np.array(ampLst)

pRef=np.array(phaseLstRef)
aRef=np.array(ampLstRef)

freq=np.linspace(-5,5,num=n)

freqFitP=np.loadtxt(os.path.join("Giovanni","phase-fit-control-off.txt"))[:,0]+10.
pfitOff=np.loadtxt(os.path.join("Giovanni","phase-fit-control-off.txt"))[:,1]
pfitOff_ref=np.loadtxt(os.path.join("Giovanni","phase-fit-control-off.txt"))[:,2]
pfitOn=np.loadtxt(os.path.join("Giovanni","phase-fit-control-on.txt"))[:,1]

freqFitO=np.loadtxt(os.path.join("Giovanni","OD-fit-control-off.txt"))[:,0]+10.+0.57
OfitOff=np.loadtxt(os.path.join("Giovanni","OD-fit-control-off.txt"))[:,1]
OfitOn=np.loadtxt(os.path.join("Giovanni","OD-fit-control-on.txt"))[:,1]

phase_offset = np.abs(pfitOff[0]-pfitOff_ref[0])


fig=plt.figure(1, figsize=(8, 9), dpi=80,)
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)

ax0 = fig.add_subplot(211)
ax0.clear()
ds=np.linspace(4,16,num=n)
dsf=np.linspace(4,16,num=len(freqFitO))
#dsf2=np.linspace(4,16,num=len(freqFitO))
ax0.plot(ds,2*np.log(aRef[:,0]/a[:,0]),'ro',label="$\Omega_c=$???")
ax0.plot(ds,2*np.log(aRef[:,1]/a[:,1]),'bo',label="$\Omega_c=$0")
#ax0.plot(dsf,OfitOff,"b")
#ax0.plot(dsf,OfitOn,"r")

ax0.plot(freqFitP,OfitOff,"b")
ax0.plot(freqFitP,OfitOn,"r")
ax0.plot(dsf,od_new(dsf*1e6*2*np.pi/Gamma_e),"g")
ax0.plot(dsf,od0_new(dsf*1e6*2*np.pi/Gamma_e),"k")

ax0.tick_params(axis="x",which="both",labelbottom="off")
ax0.set_ylabel("OD")
trans = ax0.get_xaxis_transform() # x in data untis, y in axes fraction
ax0.annotate('(a)', xy=(4.3,0.9 ), xycoords=trans)


ax1 = fig.add_subplot(212)
ax1.clear()


dsf=np.linspace(4,16,num=len(freqFitP))
ax1.plot(ds,(p[:,0]-pRef[:,0]),'ro',label="control on")
ax1.plot(ds,(p[:,1]-pRef[:,1]),'bo',label="control off")
#ax1.plot(dsf,pfitOff,"b")
#ax1.plot(dsf,pfitOn,"r")

ax1.plot(freqFitP,pfitOff,"b")
ax1.plot(freqFitP,pfitOn,"r")
ax1.plot(dsf,ph_new(dsf*1e6*2*np.pi/Gamma_e)+0*phase_offset,"g")
ax1.plot(dsf,ph0_new(dsf*1e6*2*np.pi/Gamma_e)+0*phase_offset,"k")

plt.show()

D_phi=ph_new(Delta_s)-ph0_new(Delta_s)

D_phi_blockaded = D_phi*2*R_b/L

print("Expected cond. phase shift: %.3f rad"%D_phi_blockaded)

#int_array=[]
#p_array=[]
#p2_array=[]
#t_array=[]
#t_c_array=[]
#d_cs=np.linspace(1.0,5.,40)
#for d in d_cs:
##    i,p2=cond_phase(d,Omega_c,1)
#    i,p=cond_phase(d,Omega_c,0.)
#    int_array.append(-1*i)
#    p_array.append(p)
#
#    d_s=np.abs(i)
#    Ga=Omega_c**2/(4*d_s)
#    Ga=np.true_divide(Omega_c**2*d_s,1+np.sqrt((4*d_s**2+1)**2-4*d_s**2*1))
#    print(d_s)
#    R_b=1.0*(C6/(hbar*Gamma_e)/Ga)**(1./6.)    
##    print R_b
#    chi_nb = susceptibility(i, d, gamma_21, Omega_c)
#    chi_0b = susceptibility(i, d, gamma_21, 0)
#    p2=omega_s/(2*c) * 2*R_b * chi_0 * (np.real(chi_0b)-np.real(chi_nb))
#    p2_array.append(p2)
#    
#    chi = susceptibility(i, d, gamma_21, Omega_c)
#    trans = np.exp(- omega_s/c * L * chi_0* np.imag(chi))
##    trans = omega_s/c * L * chi_0* np.imag(chi)
#    t_array.append(trans)
#    
#    _,t_c=cond_trans(d,Omega_c,0.)
#    t_c_array.append(t_c)
#    
#    detuning = np.linspace(-10.5, -0.5, 400)
#    chi = susceptibility(detuning, d, gamma_21, Omega_c)
#    chi_2L = susceptibility(detuning, d, gamma_21, 0)
#    trans = np.exp(- omega_s/c * L * chi_0* np.imag(chi))
#    trans_2L = np.exp(- omega_s/c * L * chi_0* np.imag(chi_2L))
#    
#    plt.plot(detuning,trans)
#    plt.plot(detuning,trans_2L,c='k')
#    plt.axvline(i)