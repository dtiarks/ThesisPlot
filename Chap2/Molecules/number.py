# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 15:28:20 2017

@author: daniel
"""

#import Tomography as tom
#import quPy as qp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import io
import scipy.constants as co
import rydpy

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

def classicalRadius(n,l=0,j=0.5):
    return (ec**2/(4*np.pi*epsilon_0*co.h*rydpy.IonEnergy(n,l,j,units="Hz")))
    
def orbitalRadius(n,l=0,j=0.5):
    r,rrRR,_,_=rydpy.Numerov(n,l,j)
    imax=np.argmax(rrRR)
    return r[imax]*a0

#rcl1=classicalRadius(n1)
rcl1=orbitalRadius(n1)
Vcl1=4/3.*np.pi*np.abs(rcl1)**3

rcl2=orbitalRadius(n2)
Vcl2=4/3.*np.pi*np.abs(rcl2)**3

rBend=(8*np.pi*epsilon_0*hbar**2*(n1-rydpy.QuantumDefect(n1,0,0.5))**2)/(1*ec**2*m_e)
print rBend/a0



#plt.plot(r,rrRR)
#plt.axvline(np.abs(r[imax]))


#print rho_peak*Vcl

rhos=np.linspace(0,3e12/1e-6,10)

plot_dict={}

# plot it
f=plt.figure()


#h=pd.DataFrame(index=detuning,data=np.imag(chi))
#h2=pd.DataFrame(index=detuning,data=np.imag(chi_2L))
#plot_dict['221']={
#    'A':{'type':'plot','y':h[0].to_json(),'ylabel':u'Im($\chi / \chi_0$)','margin':(0,0.0),'num':'a'},                
#    'B':{'type':'plot','y':h2[0].to_json(),'margin':(0,0.0)},
#    'C':{'type':'axv','y':0}
#}

plt.subplot(111)
plt.ylabel(u'$\Delta \phi$')
plt.plot(rhos/1e18,rhos*Vcl1)
plt.plot(rhos/1e18,rhos*Vcl2)
#plt.xscale('log')
plt.axhline(1)
#plt.yticks(np.arange(-0.4, 0.6, 0.2))
#plt.margins(0,0.1)
#plt.xlabel(u'$\varrho$ $(1e12\ cm^-3)$')
#plt.tight_layout()

h=pd.DataFrame(index=rhos/1e18,data=rhos*Vcl1)
h2=pd.DataFrame(index=rhos/1e18,data=rhos*Vcl2)
plot_dict['111']={
    'A':{'type':'plot','y':h[0].to_json(),'ylabel':r'$\left\langle N_r\right\rangle$','xlabel':r'$\varrho$ $(10^{12}\ \mathrm{cm}^{-3})$','label':'$n=69$'},                
    'B':{'type':'plot','y':h2[0].to_json(),'label':'$n=100$'},
    'C':{'type':'axh','y':1}
}




#with io.open('avg_number.json', 'w+') as f:
#  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()