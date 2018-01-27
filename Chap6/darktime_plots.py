# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 12:31:22 2018

@author: daniel
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import io
import scipy.optimize as opt
from scipy.optimize import curve_fit


def trilobite_decay(t,alpha,omega,tau,A0):
    return A0*np.exp(-t**1/(tau**1))*(np.cos(alpha)**4+np.sin(alpha)**4+2*np.cos(alpha)**2*np.sin(alpha)**2*np.cos(omega*t))
    
def alphaMu(rat):
    alpha=1/(1+(rat)**2)
    beta=1/(1+(1/rat)**2)
    return alpha,beta
    
def etaStephan(t, rat, omega, tau, A0):
    eps = omega*t/2.
    alpha, mu = alphaMu(rat)
    return A0*np.exp(-t**1/(tau**1))*(alpha**2+mu**2+2*alpha*mu*np.cos(2*eps))
    
def phiStephan(t, rat, omega, tau, A0):
    eps = omega*t/2
    alpha, mu = alphaMu(rat)
    phi=np.angle(alpha*np.exp(-1j*eps)+mu*np.exp(1j*eps))
    phi=-(np.arctan(((mu-alpha)/(mu+alpha))*np.tan(eps))-np.pi*np.sign(mu-alpha)*np.mod(eps/np.pi+0.5,1))
    return phi
    
def phitTotal(t, Dfour, rat, omega, tau, A0):
    phiMol = phiStephan(t, rat, omega, tau, A0)
    return phiMol+2*np.pi*Dfour*t


gs_effs = np.loadtxt("gs_effs.tsv")
ryd_effs = np.loadtxt("ryd_effs.tsv")

azis = np.loadtxt("azi.tsv")

    
popt1=[0.7,2*np.pi*0.2,2.,0.2]
#popt1, pcov1 = curve_fit(trilobite_decay, ryd_effs[1:,0],ryd_effs[1:,1], popt1)
popt1, pcov1 = curve_fit(etaStephan, ryd_effs[1:,0],ryd_effs[1:,1], popt1)

#popt2=[-1,0.5]
popt2=[-0.,np.pi/2]
phiFit = lambda t, Dfour, phi1: phitTotal(t, Dfour, *popt1)+np.pi/2-phi1
popt2, pcov2 = curve_fit(phiFit, azis[:,0], azis[:,1], popt2)

errs = np.sqrt(np.diag(pcov1))
errs2 = np.sqrt(np.diag(pcov2))

print "Binding energy %.2f"%(1e3*popt1[1]/(2.*np.pi))
#print 2.*np.pi/popt1[1]
print "r: %.2f"%popt1[1]
print "tau: %.2f +- %.2f"%(popt1[2],errs[2])
print "phi0: %.2f +- %.2f"%(popt2[0],errs2[0])
print "Four photon det: %.2f +- %.2f"%(popt2[1]/(2*np.pi),errs2[1]/(2*np.pi))

ts=np.linspace(ryd_effs[0,0],ryd_effs[-1,0], 200)

plot_dict={}

f=plt.figure()

plt.subplot(121)
plt.errorbar(gs_effs[:,0], gs_effs[:,1],yerr=gs_effs[:,2], ls="", marker='o',c='b')
plt.errorbar(ryd_effs[:,0], ryd_effs[:,1],yerr=ryd_effs[:,2], ls="", marker='o',c='r')
plt.plot(ts, etaStephan(ts, *popt1), ls='-',lw=1.5,c='b')

h=pd.DataFrame(index=gs_effs[1:,0],data=gs_effs[1:,1])
h_err=pd.DataFrame(index=gs_effs[1:,0],data=gs_effs[1:,2])

h2=pd.DataFrame(index=ryd_effs[1:,0],data=ryd_effs[1:,1])
h2_err=pd.DataFrame(index=ryd_effs[1:,0],data=ryd_effs[1:,2])

h3=pd.DataFrame(index=ts,data=etaStephan(ts, *popt1))

plot_dict['121']={
    'A':{'type':'errorbar','y':h[0].to_json(),'yerr':h_err[0].to_json(),'ylabel':u'Ausleseeffizienz','xlabel':u'Dunkelzeit ($\mu s$)','num':'a','xlim':(0.,8.),'ylim':(0,0.15),'label':'R'},
    'B':{'type':'errorbar','y':h2[0].to_json(),'yerr':h2_err[0].to_json(),'label':'L'},
    'C':{'type':'plot','y':h3[0].to_json()}
}


plt.subplot(122)
plt.plot(ts, phiFit(ts, *popt2), ls='-',lw=1.5,c='b')
plt.errorbar(azis[:,0], azis[:,1], yerr=azis[:,2], ls="", marker='o',lw=1.5,c='k')

h=pd.DataFrame(index=azis[1:,0],data=azis[1:,1])
h_err=pd.DataFrame(index=azis[1:,0],data=azis[1:,2])

h2=pd.DataFrame(index=ts,data=phiFit(ts, *popt2))

plot_dict['122']={
    'A':{'type':'errorbar','y':h[0].to_json(),'yerr':h_err[0].to_json(),'ylabel':u'Azimuth (rad)','xlabel':u'Dunkelzeit ($\mu s$)','num':'b','xlim':(0.,8.),'ylim':(-3,4)},
    'B':{'type':'plot','y':h2[0].to_json()}
}

with io.open('darktime.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))

plt.show()
