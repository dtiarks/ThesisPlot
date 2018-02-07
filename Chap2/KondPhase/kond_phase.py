# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 11:33:26 2017

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
from scipy.optimize import newton
from scipy import integrate

c = 299792458 # m/s, speed of light CODATA 2014
a0 = 0.52917721067e-10 # m, Bohr radius
C6 = 2.3e23 * 4.36e-18 * a0**6 # Jm^6, Van-der-Waals coefficient for the 67s - 69s
hbar = 6.626070040e-34/(2 * np.pi) # Js, Planck constant, CODATA 2014
rho_peak = 2.0e12/1e-6 # peak density in cm^-3/centi^-3
d = 2.534e-29 # Cm, dipole matrix element (D. A. Steck)
Gamma_e = 2*np.pi * 6.065e6 # decay rate (D. A. Steck)
#Gamma_e = 2*np.pi * 5.75e6 # decay rate (D. A. Steck)
epsilon_0 = 8.854187817e-12 # dielectric constant, CODATA 2014
L = 61e-6 # medium length in m
#L = 48e-6 # medium length in m
omega_s = 2*np.pi * 384.23e12 # rad/s, transition frequency
#omega_s = 2*np.pi * 377e12 # rad/s, transition frequency
gamma_21 = 0.0577314
#gamma_21 = 0.00001*1.3/5.75
chi_0 = 2*rho_peak*d**2 / (epsilon_0*hbar*Gamma_e) # prefactor of the susceptibility for the cycling transition (|R> polarization)
R_b=18e-6
#chi_0 = c*9./(omega_s*L)

def susceptibility(Delta_s, Delta_c, gamma_21, Omega_c, ladder=True):
    delta = (Delta_s + (-1 + 2*int(ladder)) * Delta_c)
    return 1j*(gamma_21 - 2j * delta)/(np.abs(Omega_c)**2 + (1 - 2j * Delta_s)*(gamma_21 - 2j * delta))

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
#    intersection=-d_c
    
    im_chi_vdw = lambda x: np.imag(susceptibility(intersection, d_c - vdW_pot(x, 1e-10)/(hbar*Gamma_e), gamma_21, om_c))
    t_1=omega_s/c *  chi_0 *integrate.quad(im_chi_vdw, -L/2, L/2)[0]
    
    
    return intersection,np.exp(-t_1)
#    return intersection,t_1
    
    
# the parameters, all in units of Gamma_3
Delta_c = 2*np.pi*15.*10**6/Gamma_e
#Delta_c = 2*np.pi*8.*10**6/Gamma_e

Omega_c = 2*np.pi*10.*10**6/Gamma_e
#Omega_c = 2*np.pi*5.6*10**6/Gamma_e

int_array=[]
p_array=[]
p2_array=[]
t_array=[]
t_c_array=[]
d_cs=np.linspace(1.0,5.,40)
for d in d_cs:
    i,p2=cond_phase(d,Omega_c,1)
#    i,p=cond_phase(d,Omega_c,0.)
    int_array.append(-1*i)
    p_array.append(p2)

    d_s=np.abs(i)
    Ga=Omega_c**2/(4*d_s)
    Ga=np.true_divide(Omega_c**2*d_s,1+np.sqrt((4*d_s**2+1)**2-4*d_s**2*1))
#    print(d_s)
    R_b=1.0*(C6/(hbar*Gamma_e)/Ga)**(1./6.)    
#    print R_b
    chi_nb = susceptibility(i, d, gamma_21, Omega_c)
    chi_0b = susceptibility(i, d, gamma_21, 0)
    p2=omega_s/(2*c) * 2*R_b * chi_0 * (np.real(chi_0b)-np.real(chi_nb))
    p2_array.append(p2)
    
    chi = susceptibility(i, d, gamma_21, Omega_c)
    trans = np.exp(- omega_s/c * L * chi_0* np.imag(chi))
#    trans = omega_s/c * L * chi_0* np.imag(chi)
    t_array.append(trans)
    
    _,t_c=cond_trans(d,Omega_c,1.)
    t_c_array.append(t_c)
    
    detuning = np.linspace(-10.5, -0.5, 400)
    chi = susceptibility(detuning, d, gamma_21, Omega_c)
    chi_2L = susceptibility(detuning, d, gamma_21, 0)
    trans = np.exp(- omega_s/c * L * chi_0* np.imag(chi))
    trans_2L = np.exp(- omega_s/c * L * chi_0* np.imag(chi_2L))
    
    plt.plot(detuning,trans)
    plt.plot(detuning,trans_2L,c='k')
    plt.axvline(i)



# plot it
f=plt.figure()


#h=pd.DataFrame(index=detuning,data=np.imag(chi))
#h2=pd.DataFrame(index=detuning,data=np.imag(chi_2L))
#plot_dict['221']={
#    'A':{'type':'plot','y':h[0].to_json(),'ylabel':u'Im($\chi / \chi_0$)','margin':(0,0.0),'num':'a'},                
#    'B':{'type':'plot','y':h2[0].to_json(),'margin':(0,0.0)},
#    'C':{'type':'axv','y':0}
#}

int_array=-1.*np.array(int_array)

plt.subplot(121)
plt.plot(int_array , p_array, ls='-',lw=1.5,c='r')
plt.plot(int_array , p2_array, ls='-',lw=1.5,c='b')
plt.axhline(np.pi)
#plt.xlabel(u'$\Delta_s / \Gamma_3$')
plt.ylabel(u'$\Delta \phi$')
#plt.yticks(np.arange(-0.4, 0.6, 0.2))
#plt.margins(0,0.1)
plt.xlabel(u'$\Delta_s$ $(\Gamma_e)$')
#plt.tight_layout()

plot_dict={}

h=pd.DataFrame(index=int_array,data=p_array)
h2=pd.DataFrame(index=int_array,data=p2_array)
plot_dict['121']={
    'A':{'type':'plot','y':h[0].to_json(),'ylabel':u'$\Delta \\varphi$ (rad)','xlabel':u'$\Delta_s/ \Gamma_3$','margin':(0.1,0.2),'num':'a','xlim':(-4.2,-1.6),'ylim':(1,6),'label':'Numerisch'},                
    'B':{'type':'plot','y':h2[0].to_json(),'label':'N\"aherung'},
    'C':{'type':'axh','y':np.pi}
}

# detuned eit
# the parameters, all in units of Gamma_3
#Delta_c = 1.5
#Delta_p = -1.6
#gamma_21 = 0
#Omega_c = 1.5


h2=pd.DataFrame(index=int_array,data=t_array)
h=pd.DataFrame(index=int_array,data=t_c_array)
plot_dict['122']={
    'A':{'type':'plot','y':h[0].to_json(),'xlabel':u'$\Delta_s/ \Gamma_3$','ylim':(0.00,0.9),'xlim':(-4.2,-1.6),'num':'b','ylabel':u'Transmission'},
    'B':{'type':'plot','y':h2[0].to_json()},
    'C':{'type':'axh','y':-1}
}


plt.subplot(122)
plt.plot(int_array , t_array, ls='-',lw=1.5,c='r')
plt.plot(int_array , t_c_array, ls='-',lw=1.5,c='b')
#plt.plot(Rs*1e6 , np.real(chi), ls='-',lw=1.5)
#plt.xlabel(u'R $(\mu m)$')
#plt.ylabel(u'Im($\chi / \chi_0$)')
#plt.margins(0,0.1)
#plt.axvline(-1.5, ls='--',lw=1.0)
#plt.tight_layout()



with io.open('cond_phase.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))
#
#
plt.show()
#
#i,p=cond_phase(Delta_c,Omega_c,0.)
#
#print("Num: %.3f"%p)
#
#
#d=Delta_c
#d_s=Delta_c
#
#Ga=Omega_c**2/(4*d_s)
#Ga=np.true_divide(Omega_c**2*d_s,1+np.sqrt((4*d_s**2+1)**2-4*d_s**2*1))
#R_b=1.0*(C6/(hbar*Gamma_e)/Ga)**(1./6.)    
#
#chi_nb = susceptibility(i+2*np.pi*0.57*10**6/Gamma_e, d, gamma_21, Omega_c)
#chi_0b = susceptibility(i, d, gamma_21, 0)
#p2=omega_s/(2*c) * 2*R_b * chi_0 * (np.real(chi_0b)-np.real(chi_nb))
#
#print("Approx: %.3f"%p2)
#
#plt.show()
#plot_dict={}
