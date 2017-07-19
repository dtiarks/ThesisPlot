# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 11:29:40 2017

@author: daniel
"""

import Tomography as tom
import quPy as qp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import io
import scipy.constants as co
#import rydpy

c = 299792458 # m/s, speed of light CODATA 2014
a0 = 0.52917721067e-10 # m, Bohr radius
C6 = 2.3e23 * 4.36e-18 * a0**6 # Jm^6, Van-der-Waals coefficient for the 67s - 69s
hbar = 6.626070040e-34/(2 * np.pi) # Js, Planck constant, CODATA 2014
rho_peak = 2.0e12/1e-6 # peak density in cm^-3/centi^-3
d = 2.534e-29 # Cm, dipole matrix element (D. A. Steck)
Gamma_e = 2*np.pi * 6.065e6 # decay rate (D. A. Steck)
epsilon_0 = 8.854187817e-12 # dielectric constant, CODATA 2014
L = 61e-6 # medium length in m
omega_s = 2*np.pi * 384.23e12 # rad/s, transition frequency
gamma_21 = 0.0577314
chi_0 = 2*rho_peak*d**2 / (epsilon_0*hbar*Gamma_e) # prefactor of the susceptibility for the cycling transition (|R> polarization)
R_b=18e-6
n1=69
n2=100
ec=co.e
m_e=co.electron_mass

#def classicalRadius(n,l=0,j=0.5):
#    return (ec**2/(4*np.pi*epsilon_0*co.h*rydpy.IonEnergy(n,l,j,units="Hz")))
#    
#def orbitalRadius(n,l=0,j=0.5):
#    r,rrRR,_,_=rydpy.Numerov(n,l,j)
#    imax=np.argmax(rrRR)
#    return r[imax]*a0

def moleculeState(t,pm,omega):
    return np.cos(pm)*np.array([1,0])+np.sin(pm)*np.exp(-1j*omega*t)*np.array([0,1])
    
def nuMolecule(t,pm,omega):
    return np.cos(pm)**4 + np.sin(pm)**4 + 2*np.sin(pm)**2*np.cos(pm)**2*np.cos(omega*t)
    
def nuNumerical(ts,pm,omega):
    s=np.array([np.inner(moleculeState(0,pm,omega),moleculeState(tx,pm,omega)) for tx in ts])
    return np.abs(s)**2
    
def phiNumerical(ts,pm,omega):
    s=np.array([np.angle(moleculeState(tx,pm,omega)[0]+moleculeState(tx,pm,omega)[1]) for tx in ts])
    return s
    
#rcl1=classicalRadius(n1)
#rcl1=orbitalRadius(n1)
#Vcl1=4/3.*np.pi*np.abs(rcl1)**3
#
#rcl2=orbitalRadius(n2)
#Vcl2=4/3.*np.pi*np.abs(rcl2)**3
#
#rBend=(8*np.pi*epsilon_0*hbar**2*(40-rydpy.QuantumDefect(40,0,0.5))**2)/(1*ec**2*m_e)
#print rBend/a0

pM1=np.pi/12.
pM2=np.pi/6.
omega=2*np.pi*220e3
ts=np.linspace(0,10e-6,100)
#nu=np.abs(np.outer(moleculeState(np.zeros_like(ts),pM,0),moleculeState(ts,pM,omega)))**2

plot_dict={}

# plot it
f=plt.figure()


h=pd.DataFrame(index=1e6*ts,data=nuNumerical(ts,pM1,omega))
h2=pd.DataFrame(index=1e6*ts,data=nuNumerical(ts,pM2,omega))
plot_dict['121']={
    'A':{'type':'plot','y':h[0].to_json(),'ylabel':u'$\eta_L(t)$','xlabel':r'Speicherzeit $(\mu s)$','ylim':(0,1.2),'num':'a','label':r'$\alpha=\pi/12$'},                
    'B':{'type':'plot','y':h2[0].to_json(),'label':r'$\alpha=\pi/6$'}
}

plt.subplot(121)
plt.ylabel(u'$\Delta \phi$')
plt.plot(1e6*ts,nuNumerical(ts,pM1,omega))
plt.plot(1e6*ts,nuNumerical(ts,pM2,omega))
#plt.plot(1e6*ts,nuMolecule(ts,pM,omega))
#plt.axhline(1)

h=pd.DataFrame(index=1e6*ts,data=phiNumerical(ts,pM1,omega))
h2=pd.DataFrame(index=1e6*ts,data=phiNumerical(ts,pM2,omega))
plot_dict['122']={
    'A':{'type':'plot','y':h[0].to_json(),'ylabel':r'$\Delta\varphi_{\mathrm{mol}}$','xlabel':r'Speicherzeit $(\mu s)$','label':r'$\alpha=\pi/12$','num':'b','ylim':(-1.85,0.8)},                
    'B':{'type':'plot','y':h2[0].to_json(),'label':r'$\alpha=\pi/6$'}
}

plt.subplot(122)
plt.ylabel(u'$\Delta \phi$')
plt.plot(1e6*ts,phiNumerical(ts,pM1,omega))
plt.plot(1e6*ts,phiNumerical(ts,pM2,omega))


with io.open('memory.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()