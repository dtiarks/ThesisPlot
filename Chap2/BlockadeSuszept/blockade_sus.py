# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 16:51:51 2017

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

c = 299792458 # m/s, speed of light CODATA 2014
a0 = 0.52917721067e-10 # m, Bohr radius
C6 = 2.3e23 * 4.36e-18 * a0**6 # Jm^6, Van-der-Waals coefficient for the 67s - 69s
hbar = 6.626070040e-34/(2 * np.pi) # Js, Planck constant, CODATA 2014
rho_peak = 2.0e12/1e-6 # peak density in cm^-3/centi^-3
d = 2.534e-29 # Cm, dipole matrix element (D. A. Steck)
Gamma_e = 2*np.pi * 6.065e6 # decay rate (D. A. Steck)

def susceptibility(Delta_s, Delta_c, gamma_21, Omega_c, ladder=True):
    delta = (Delta_s + (-1 + 2*int(ladder)) * Delta_c)
    return 1j*(gamma_21 - 2j * delta)/(np.abs(Omega_c)**2 + (1 - 2j * Delta_s)*(gamma_21 - 2j * delta))

def vdW_pot(r, r0):
    return -C6 * (r-r0)**-6

# the parameters, all in units of Gamma_3
Delta_c = 0
Delta_p = 0
gamma_21 = 0
Omega_c = 1.8

# calculate the curve
Rs = np.linspace(0, 40e-6, 400)
chi = susceptibility(Delta_p, Delta_c-vdW_pot(Rs, 1e-10)/(hbar*Gamma_e), gamma_21, Omega_c)
#chi_2L = susceptibility(detuning, Delta_c, gamma_21, 0)

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

plt.subplot(121)
plt.plot(Rs*1e6 , np.real(chi), ls='-',lw=1.5)
plt.plot(Rs*1e6 , np.imag(chi), ls='-',lw=1.5)
#plt.xlabel(u'$\Delta_s / \Gamma_3$')
plt.ylabel(u'$\chi / \chi_0$')
#plt.yticks(np.arange(-0.4, 0.6, 0.2))
#plt.margins(0,0.1)
plt.axvline(0, ls='--',lw=1.0)
plt.xlabel(u'R $(\mu m)$')
#plt.tight_layout()

h=pd.DataFrame(index=Rs*1e6,data=np.real(chi))
h2=pd.DataFrame(index=Rs*1e6,data=np.imag(chi))
plot_dict['121']={
    'A':{'type':'plot','y':h[0].to_json(),'ylabel':u'$\chi / \chi_0$','xlabel':u'R $(\mu m)$','margin':(0,0.0),'num':'a','label':u'Re($\chi$)'},                
    'B':{'type':'plot','y':h2[0].to_json(),'margin':(0,0.025),'label':u'Im($\chi$)'}
#    'C':{'type':'axv','y':0}
}

# detuned eit
# the parameters, all in units of Gamma_3
Delta_c = 1.5
Delta_p = -1.6
gamma_21 = 0
Omega_c = 1.5

# calculate the curve
chi = susceptibility(Delta_p, Delta_c-vdW_pot(Rs, 1e-10)/(hbar*Gamma_e), gamma_21, Omega_c)
chi_2L = susceptibility(Delta_p, Delta_c, gamma_21, 0)
chi_EIT = susceptibility(Delta_p, Delta_c, gamma_21, Omega_c)

h=pd.DataFrame(index=Rs*1e6,data=np.real(chi))
h2=pd.DataFrame(index=Rs*1e6,data=np.imag(chi))
plot_dict['122']={
    'A':{'type':'plot','y':h[0].to_json(),'xlabel':u'R $(\mu m)$','margin':(0,0.0),'num':'b','label':u'Re($\chi$)'},                
    'B':{'type':'plot','y':h2[0].to_json(),'margin':(0,0.025),'label':u'Im($\chi$)'}
#    'C':{'type':'axv','y':0}
}


plt.subplot(122)
plt.plot(Rs*1e6 , np.real(chi), ls='-',lw=1.5)
plt.plot(Rs*1e6 , np.imag(chi), ls='-',lw=1.5)
plt.axhline(np.imag(chi_2L))
plt.axhline(np.imag(chi_EIT))
plt.xlabel(u'R $(\mu m)$')
#plt.ylabel(u'Im($\chi / \chi_0$)')
#plt.margins(0,0.1)
#plt.axvline(-1.5, ls='--',lw=1.0)
#plt.tight_layout()



with io.open('blockade_suszept.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()