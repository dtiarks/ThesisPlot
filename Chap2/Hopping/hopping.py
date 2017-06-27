# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 16:34:55 2017

@author: daniel
"""

import matplotlib.pyplot as plt
from scipy.optimize import newton
import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import newton
import numpy as np
import datetime

from scipy.integrate import odeint
from scipy.misc import derivative


### PARAMETERS AND CONSTANTS

hbar = 6.626070040e-34/(2 * np.pi) # Js, Planck constant, CODATA 2014
rho_peak = 2.0e12/1e-6 # peak density in cm^-3/centi^-3
d = 2.534e-29 # Cm, dipole matrix element (D. A. Steck)
Gamma_e = 2*np.pi * 6.065e6 # decay rate (D. A. Steck)
epsilon_0 = 8.854187817e-12 # dielectric constant, CODATA 2014
L = 60e-6 # medium length in m
omega_s = 2*np.pi * 384.23e12 # rad/s, transition frequency
c = 299792458 # m/s, speed of light CODATA 2014
a0 = 0.52917721067e-10 # m, Bohr radius
C6 = 2.3e23 * 4.36e-18 * a0**6 # Jm^6, Van-der-Waals coefficient for the 67s - 69s
# interaction potential V=-C6/R^6, converted from atomic units with the Hartree 
# energy and the Bohr radius

# the parameters, all in units of Gamma_3
Delta_c = 2*np.pi*15.0*10**6/Gamma_e
gamma_21 = 0.0577314
Omega_c = 2*np.pi*8.*10**6/Gamma_e

def susceptibility(Delta_s, Delta_c, gamma_21, Omega_c, ladder=True):
    delta = (Delta_s + (-1 + 2*int(ladder)) * Delta_c) # two photon detuning
    return 1j*(gamma_21 - 2j * delta)/(np.abs(Omega_c)**2 + (1 - 2j * Delta_s)*(gamma_21 - 2j * delta))

chi_0 = 2*rho_peak*d**2 / (epsilon_0*hbar*Gamma_e) # prefactor of the susceptibility for the cycling transition (|R> polarization)
intersection = newton(lambda x: np.imag(susceptibility(x, Delta_c, gamma_21, Omega_c) - susceptibility(x, Delta_c, gamma_21, 0)), -Delta_c)



def vdW_pot(r, r0):
    return -C6 * (r-r0)**-6

def index(det):
    chi=chi_0 * susceptibility(intersection, det, gamma_21, Omega_c)
    n=np.sqrt(1+np.real(chi)) # index of refraction
    return n

def group_velocity(d_c):
    ''' group velocities in meters/second for the given susceptibility array chi
    with frequency distance d_omega between the points.'''
    d_o = Delta_c*0.01
    n=index(d_c)
    dn=derivative(index, d_c, dx=d_o)
    v_gr =  c/(n + omega_s * dn)
    return v_gr# formula from Fleischhauer's Rev. Mod. Phys. 2005  

# calculate the intersection of the imaginary parts


#print "differential phase shift over whole cloud at equal absorption:"
#print omega_s/(2*c) * L * chi_0 *np.real(susceptibility(intersection, Delta_c, gamma_21, Omega_c)- susceptibility(intersection, Delta_c, gamma_21, 0))


# calculate the transmission and phase curves
detuning = np.linspace(-4.5, -0.5, 400)

R0=-L/2
t=np.linspace(0,0.63e-6*Gamma_e,50000)

def func(R,ts):
    d_c = Delta_c - vdW_pot(R, 0)/(hbar*Gamma_e)
    d_o = (detuning[1]-detuning[0]) * Gamma_e
    d_o = Delta_c*0.01
    n=index(d_c)
    dn=derivative(index, d_c, dx=d_o)
    v_gr =  c/(n + omega_s * dn)
    return v_gr
    
Rs=odeint(func,R0,t)

D_Ep=5e6

def PRR(tss,Rss):
#    v=np.array(vdW_pot(Rss,0)/(hbar))[:,0]
    v=np.array(vdW_pot(Rss,0)/(hbar))
    A=np.array(v**2/(v**2+D_Ep**2))
    Om=2*np.pi*np.sqrt(v**2+D_Ep**2)
    P=A*np.sin(0.5*Om*tss)**2
    return P

fig=plt.figure(0,figsize=(9,10))

plt.suptitle("Hopping dynamics",fontsize=15, fontweight='bold')

ax1=fig.add_subplot(211)
ax1.plot(1e6*t/Gamma_e,1e6*Rs,label='w/ excitation')
ax1.plot(1e6*t/Gamma_e,1e6*group_velocity(Delta_c)*t+1e6*R0,label='w/o excitation')
ax1.axhline(-1e6*R0,c='k',label='Medium border')
ax1.axhline(1e6*R0,c='k')
ax1.axhline(-15,c='r',label=u'$R_b=15\mu m$')
ax1.axhline(15,c='r')
ax1.set_ylabel("Distance (um)")
ax1.set_xlabel("Time (us)")
ax1.legend(loc=2)

ax2=fig.add_subplot(212)
#ax2.plot(1e6*t/Gamma_e,PRR(t/Gamma_e,Rs))
RbMask=np.abs(Rs)>15e-6
Rbs=Rs[RbMask]
ax2.plot(1e6*Rbs,PRR(np.array(t)[RbMask[:,0]]/Gamma_e,Rbs))
#ax2.axvline(0.0814,c='r')
#ax2.axvline(-0.083,c='r')
ax2.set_ylabel("Population reversed state")
ax2.set_xlabel("Distance travalled (um)")

plt.savefig("hopping.pdf")
plt.show()



