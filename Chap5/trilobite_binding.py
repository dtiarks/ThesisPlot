# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 12:27:15 2017

@author: daniel
"""

import numpy as np
import scipy.constants as const
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import rydpy2 as rydpy

Ry = 109736.605e2# * const.c #Rydberg constant (Hz) from [1]
a0 = 0.52917721067e-10 # m, Bohr radius
hbar = 6.626070040e-34/(2 * np.pi) # Js, Planck constant, CODATA 2014
Eh=const.physical_constants["Hartree energy"][0]
me=const.m_e
u=const.physical_constants["atomic mass constant"][0]
m=43.5*u



bound=8910
bound_u=150
#for n in np.arange(69,70):
r, rrRR, Vscat, G = rydpy.Numerov(69,0.,1./2)
rV,Vs2,_=rydpy.ScatteringPotential(69,0.,1./2)

V_min_r=np.argmin(Vs2[bound:])
V_min=np.min(Vs2[bound:])

V_frac=0.7

VscatMin=Vs2[V_min_r+bound]
iStart=V_min_r+bound-1
while np.abs(Vs2[iStart]) > np.abs(VscatMin)*V_frac:
    iStart=iStart-1

print iStart
iStop=V_min_r+bound+1
while np.abs(Vs2[iStop]) > np.abs(VscatMin)*V_frac:
    iStop=iStop+1

print iStop

off=V_min*Eh
x0=rV[bound:][V_min_r]*a0

def harmFunc(x, w):
  return 1./2*m*w**2*(x-x0)**2+off



popt1=[2*np.pi*400e3]
popt1, pcov1 = curve_fit(harmFunc, rV[iStart:iStop]*a0,Vs2[iStart:iStop]*Eh, popt1)
b_state = harmFunc(rV[iStart:iStop]*a0,*popt1)

v0=popt1[0]/2+off/hbar

print "tau: %f"%(-1*1e6*2*np.pi/(v0))
print "omega: %f"%(v0/(2.*np.pi))

fig=plt.figure(0)
fig.clear()
ax0=fig.add_subplot(211)
ax1=fig.add_subplot(212)

start=2000
ax0.plot(r[start:], rrRR[start:]/(r[start:]*r[start:]))
ax0.set_ylabel("$|r R_{n,l}|^2$")
ax0.set_xlabel("$r/a_0$")

#ax0.set_title("$%ds_{1/2}$" %n)

ax1.plot(rV[5000:]*a0, (Vs2[5000:]*Eh/hbar/(2*np.pi))/1e3)
ax1.plot(rV[iStart:iStop]*a0, (b_state/hbar/(2*np.pi))/1e3)
ax1.axvline(rV[bound:][V_min_r]*a0)
ax1.axhline((V_min*Eh/hbar/(2*np.pi))/1e3)
#  ax1.set_ylim([-5.2,0])
ax1.grid()
#  ax1.plot(p[0], harmFunc(*p)*2*Ry/1e6, label=fitStr)
ax1.legend(loc=4, prop={"size" : 9})
  #ax1.set_xlim([500,2500])
ax1.set_ylabel("$V_{scatt} (MHz)$")
ax1.set_xlabel("$r/a_0$")
fig.show()

  #fig.set_size_inches((7,6))
  #plt.savefig("example.pdf")

  #f=open("theory.tsv", "a")
  #f.write("%d\t%e\n" % (n, nu0))
  #f.close()