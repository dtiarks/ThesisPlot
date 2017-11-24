# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 15:57:42 2017

@author: daniel
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import io
import scipy.optimize as opt

freq = [49.8933, 57.7919, 55.29, 46.71, 58.315]
Temp = [25.0, 30.0, 35.0, 40.0, 31.74]
Ts =  np.arange(25.0, 40.0, 0.01)

def fitfunc(T, T0, A, B):
    return A*(T - T0)**2 + B
    
def decay(t, t0, tau, A, B):
    return A*np.exp(-(t-t0)/tau)+B
    


plot_dict={}

f=plt.figure()

spectra_dir="spectra"

trans_fit_off = np.loadtxt(os.path.join(spectra_dir, "Trans-And-Phase-Off-fit.txt"))
trans_fit_on = np.loadtxt(os.path.join(spectra_dir, "Trans-And-Phase-On-fit.txt"))

trans = np.loadtxt(os.path.join(spectra_dir, "RawSpectraWithError.txt"))

phase_fit_off = np.loadtxt(os.path.join(spectra_dir, "Wrapped-Phase-Off-With-Error.txt"))
phase_fit_on = np.loadtxt(os.path.join(spectra_dir, "Wrapped-Phase-On-With-Error.txt"))


h=pd.DataFrame(index=trans_fit_off[:,0],data=trans_fit_off[:,1])
h2=pd.DataFrame(index=trans_fit_on[:,0],data=trans_fit_on[:,1])

h3=pd.DataFrame(index=trans[:,0],data=trans[:,1])
h3_err=pd.DataFrame(index=trans[:,0],data=trans[:,2])
h4=pd.DataFrame(index=trans[:,5],data=trans[:,6])
h4_err=pd.DataFrame(index=trans[:,5],data=trans[:,7])

plot_dict['121']={
    'A':{'type':'plot','y':h[0].to_json()},
    'B':{'type':'plot','y':h2[0].to_json()},
    'C':{'type':'errorbar','y':h3[0].to_json(),'yerr':h3_err[0].to_json(),'ylabel':u'Transmission','xlabel':u'Signalverstimmung $\Delta_s/2\pi$ (MHz)','num':'a','xlim':(-25.,5.),'ylim':(0.,1.0),'label':'Daten'},
    'D':{'type':'errorbar','y':h4[0].to_json(),'yerr':h4_err[0].to_json()},
    'E':{'type':'axv','y':-10.}
    
}

plt.subplot(121)
plt.plot(trans_fit_off[:,0], trans_fit_off[:,1], ls='-',lw=1.5,c='r')
plt.plot(trans_fit_on[:,0], trans_fit_on[:,1], ls='-',lw=1.5,c='b')
plt.errorbar(trans[:,0], trans[:,1], yerr=trans[:,2], ls="", marker='o',lw=1.5,c='r')
plt.errorbar(trans[:,5], trans[:,6], yerr=trans[:,7], ls="", marker='o',lw=1.5,c='b')


h=pd.DataFrame(index=trans_fit_off[:,0],data= trans_fit_off[:,2])
h2=pd.DataFrame(index=trans_fit_on[:,0],data=trans_fit_on[:,2])

h3=pd.DataFrame(index=phase_fit_off[:,0],data=phase_fit_off[:,1])
h3_err=pd.DataFrame(index=phase_fit_off[:,0],data=phase_fit_off[:,2])
h4=pd.DataFrame(index=phase_fit_on[:,0],data=phase_fit_on[:,1])
h4_err=pd.DataFrame(index=phase_fit_on[:,0],data=phase_fit_on[:,2])

plot_dict['122']={
    'A':{'type':'plot','y':h[0].to_json()},
    'B':{'type':'plot','y':h2[0].to_json()},
    'C':{'type':'errorbar','y':h3[0].to_json(),'yerr':h3_err[0].to_json(),'ylabel':u'Phase (rad)','xlabel':u'Signalverstimmung $\Delta_s/2\pi$ (MHz)','num':'b','xlim':(-25.,5.),'ylim':(-8.,8.),'label':'Daten'},
    'D':{'type':'errorbar','y':h4[0].to_json(),'yerr':h4_err[0].to_json()},
    'E':{'type':'axv','y':-10.}
    
}

plt.subplot(122)
plt.plot(trans_fit_off[:,0], trans_fit_off[:,2], ls='-',lw=1.5,c='r')
plt.plot(trans_fit_on[:,0], trans_fit_on[:,2], ls='-',lw=1.5,c='b')
plt.errorbar(phase_fit_off[:,0], phase_fit_off[:,1], yerr=phase_fit_off[:,2], ls="", marker='o',lw=1.5,c='r')
plt.errorbar(phase_fit_on[:,0], phase_fit_on[:,1], yerr=phase_fit_on[:,2],  ls="",marker='o',lw=1.5,c='b')
plt.ylim((-8,8))



with io.open('pol_spectra.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()