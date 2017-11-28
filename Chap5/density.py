# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 15:51:13 2017

@author: daniel
"""

import numpy as np
from subprocess import call
import matplotlib.pyplot as plt
import scipy as sc
import sys
import os
from datetime import datetime
from scipy import stats
import numpy.ma as ma
from scipy.optimize import curve_fit
import scipy.constants as c
import scipy.integrate as integ

u=c.physical_constants["atomic mass constant"][0]
m=87*u
a0=c.physical_constants["Bohr radius"][0]

pxs=3.28*10**-6
pxAndor=4
TOF=12*10**-3
DF=46.6
PlugDF=2.5
om_r=2*np.pi*90
L=PlugDF*DF*10**-6
#L=120e-6
P=0.198*0.5
#P=0.1
#w1=25*10**-6
w=25*10**-6
Dp=4.41206913e-05

au=4*np.pi*c.epsilon_0*a0**3
alpha=-254.*au

I0=2*P/(np.pi*w**2)
U0=I0*(-0.5*alpha/(c.epsilon_0*c.c))

#U0=21e-6*c.k

def Iplug(z,z0):
    return I0*np.exp(-2*(z-z0)**2/w**2)

def Uplug(z):
    #return Iplug(z,z0)*(-0.5*alpha/(c.epsilon_0*c.c))
    return U0*(np.exp(-2*(z+L/2)**2/w**2)+np.exp(-2*(z-L/2)**2/w**2))

def Pre(z,T):
    return np.exp(-(Uplug(z))/(c.k*T))

def dens(z,T):
   

    return (np.exp(-((z)/Dp)**4))

def makeDensity(atomnumber,temp):
    index=range(len(atomnumber))
    rho=[]
    for i in index:
        a=1.e5
        T=524.*10**-9
        #r=0.31831*a*m*om_r**2/(c.k*T)*integ.quad(Pre,-L/2,L/2,args=(T,))
        y,err=integ.quad(Pre,-L/2,L/2,args=(T))
        #y,err=integ.quad(dens,-L,L,args=(T))
        #r=(0.31831*a*m*om_r**2/(c.k*T))/y
        r=y*2*np.pi*c.k*T/(m*om_r**2)
        
        #r=1.8128*Dp*2*np.pi*c.k*T/(m*om_r**2)
        rho.append(a/r)
	print np.sqrt(c.k*T/m)/om_r
    return np.asarray(rho)
    
def makeDensityS(atomnumber,temp):

    
        #r=0.31831*a*m*om_r**2/(c.k*T)*integ.quad(Pre,-L/2,L/2,args=(T,))
    y,err=integ.quad(Pre,-L/2,L/2,args=(temp))
    #y,err=integ.quad(dens,-L,L,args=(temp))
        #r=(0.31831*a*m*om_r**2/(c.k*T))/y
    r=y*2*np.pi*c.k*temp/(m*om_r**2)
        
        #r=1.8128*Dp*2*np.pi*c.k*T/(m*om_r**2)
    return atomnumber/r

def main(argv):
  
  
#  atomNumerbs=np.loadtxt("img_001/fit_results.had")
#  gaussFits=np.loadtxt("img_001/fit_results2.had")

  
  
  def second(at):
      return 1.31*at+1.2e+04
      
  def secondNoOffset(at):
      return 1.43*at
      
  def secondTemp(T):
      return 1.242*T-1.1e-7
      
  def secondTempNoOffset(T):
      return 1.03*T
  
#  TMean=[]
#  AtomMean=[]
#  rhoMean=[]
#  rhoErr=[]
#  
#  txMean=np.mean(gaussFits[:,3])
#  print txMean
#  #T=(m/c.k)*np.asarray(pxs*t[:,1])**2/(TOF**2)
#  T=(m/c.k)*np.asarray(pxs*gaussFits[:,3])**2/(TOF**2)
#  #print np.mean(T)
#  TMean.append(np.mean(T))
#  #TMean.append(secondTemp(T))
#  
#  Ats=secondNoOffset(  atomNumerbs[:,1])
#  print Ats
#  #AtomMean.append(At)
#  AtomMean.append(Ats)
  
  #rho=makeDensity(t[:,0],T)
  T=0.5e-6
  T_second=secondTemp(T)
  
  At=7.4e4
  At_second=secondNoOffset(At)
  rho=makeDensityS(At_second,T_second)
  
  print "Atomnumber: "
  print np.mean(At_second)
  print "Temperature: "
#  Tm=np.mean(secondTempNoOffset(T))
  print T_second
  
#  rhoMean.append(np.mean(rho))
#  rhoErr.append(stats.sem(rho))
  print "Density (cm -3): "
#  print np.asarray(rhoMean)/1e18
  print rho/1e18
#  print "Density  (cm -3): "
#  print np.asarray(rhoErr)/1e18
  
  print "Barrier height (kb):"
  print U0/c.k
  
  print "Plug beam distance:"
  print L
   
#  zs=np.linspace(-L/2,L/2,1000)
#  P=Pre(zs,Tm)
#  Ph=P[:len(P)/2]
#  fwhmind=2*np.argwhere(np.abs(Ph-0.5)<0.005)[0]
#  print "FWHM"
#  print zs[fwhmind]
#  plt.plot(zs,P)
#  plt.show()


if __name__ == "__main__":
  main(sys.argv)
    