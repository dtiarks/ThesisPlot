# -*- coding: utf-8 -*-
"""
Created on Tue Jan 03 17:44:08 2017

@author: daniel
"""

import numpy as np
import scipy.constants as c
import scipy.integrate as integ
import pylab as plt
import pandas as pd
import os
import json
import io
import scipy.constants as co

u=c.physical_constants["atomic mass constant"][0]
m=87*u
a0=c.physical_constants["Bohr radius"][0]

pxs=3.28*10**-6
pxAndor=4
TOF=12*10**-3
om_r=2*np.pi*89
dF=3.
L=42.45*dF*10**-6
L_medium=41.55*dF-27.93
P=0.24
w=14.5*10**-6
Dp=4.41206913e-05

T=700*10**-9
N_Atoms=0.8*10**5

au=4*np.pi*c.epsilon_0*a0**3
alpha=-254.*au

I0=P/(np.pi*w**2)
U0=I0*(-0.5*alpha/(c.epsilon_0*c.c))

def Iplug(z,z0):
    return I0*np.exp(-2*(z-z0)**2/w**2)

def Uplug(z):
    return U0*(np.exp(-2*(z+L/2)**2/w**2)+np.exp(-2*(z-L/2)**2/w**2))

def Pre(z,T):
    return np.exp(-(Uplug(z))/(c.k*T))
    
def secondDetectAtomNumber(at):
    return 1.31*at+1.2e+04

def secondDetectTemperature(T):
    return 1.242*T-1.1e-7
  
def getTemperatureFromFit(sigx,t_TOF,T_cal=lambda x: x):
    T=T_cal((m/c.k)*np.asarray(pxs*sigx)**2/(t_TOF**2))
    return T
    
       
def getDensity(atomnumber,temp,trap='',atom_cal=lambda x: x):
    y,err=integ.quad(Pre,-L/2,L/2,args=(temp))
    
    r=y*2*np.pi*c.k*temp/(m*om_r**2)
        
    return atom_cal(atomnumber)/r

#tMean=np.mean(fits[:,3])
#ats=np.mean(atomnumbers[:,1])
  
zs=np.linspace(-100*10**-6,100*10**-6,500)


plot_dict={}

# plot it
f=plt.figure()


h=pd.DataFrame(index=1e6*zs,data=Uplug(zs)/U0)
h2=pd.DataFrame(index=1e6*zs[100:400],data=Pre(zs[100:400],T)/Pre(0,T))
plot_dict['111']={
    'A':{'type':'plot','y':h[0].to_json(),'xlabel':r'Ort $(\mu m)$','ylim':(0,1.1),'label':r'$U_p(z)/U_0$'},                
    'B':{'type':'plot','y':h2[0].to_json(),'label':r'Dichte $n_p(z)/n0$'}
}

plt.subplot(111)
plt.plot(zs,Uplug(zs)/U0)
plt.plot(zs[100:400],Pre(zs[100:400],T)/Pre(0,T))

with io.open('density_profile.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))

plt.show()
    
