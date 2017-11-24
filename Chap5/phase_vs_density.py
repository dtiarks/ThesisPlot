# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 18:08:54 2017

@author: daniel
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import io
import scipy.optimize as opt


def fitfunc(x, A):
    return A*x
    
    


plot_dict={}

f=plt.figure()

density_dir="density"

phases = np.loadtxt(os.path.join(density_dir, "PhaseVsDensity.txt"))
densitys = np.loadtxt(os.path.join(density_dir, "Densities-25um-Waist.txt"))




h3=pd.DataFrame(index=densitys/1e18,data=phases[:,1])
h3_err=pd.DataFrame(index=densitys/1e18,data=phases[:,2])

h4=pd.DataFrame(index=densitys/1e18,data=phases[:,3])
h4_err=pd.DataFrame(index=densitys/1e18,data=phases[:,4])


popt, pcov = opt.curve_fit(fitfunc, densitys/1e18, phases[:,1])
popt2, pcov2 = opt.curve_fit(fitfunc, densitys/1e18, phases[:,3])
popt3, pcov3 = opt.curve_fit(fitfunc, densitys/1e18, phases[:,5])

ds=np.linspace(0,2.2,100)

h=pd.DataFrame(index=ds,data=fitfunc(ds,*popt))
h2=pd.DataFrame(index=ds,data=fitfunc(ds,*popt2))


plot_dict['121']={
    'A':{'type':'plot','y':h[0].to_json()},
    'B':{'type':'plot','y':h2[0].to_json()},
    'C':{'type':'errorbar','y':h3[0].to_json(),'yerr':h3_err[0].to_json(),'ylabel':u'Phase (rad)','xlabel':u'Dichte ($10^{12}\ cm^{-3}$)','num':'a','xlim':(0.,2.5),'ylim':(-5.5,1),'label':'Daten'},
    'D':{'type':'errorbar','y':h4[0].to_json(),'yerr':h4_err[0].to_json()}
}

plt.subplot(121)
plt.plot(ds, fitfunc(ds,*popt), ls='-',lw=1.5,c='r')
plt.plot(ds, fitfunc(ds,*popt2), ls='-',lw=1.5,c='b')
plt.errorbar(densitys/1e18, phases[:,1], yerr=phases[:,2], ls="", marker='o',lw=1.5,c='r')
plt.errorbar(densitys/1e18, phases[:,3], yerr=phases[:,4], ls="", marker='o',lw=1.5,c='b')


h=pd.DataFrame(index=densitys/1e18,data=phases[:,5])
h_err=pd.DataFrame(index=densitys/1e18,data=phases[:,6])
h2=pd.DataFrame(index=ds,data=fitfunc(ds,*popt3))

plot_dict['122']={
    'A':{'type':'plot','y':h2[0].to_json()},
    'B':{'type':'errorbar','y':h[0].to_json(),'yerr':h_err[0].to_json(),'xlabel':u'Dichte ($10^{12}\ cm^{-3}$)','num':'b','xlim':(0.,2.5),'ylim':(0,4.0),'label':'Daten'},
}

plt.subplot(122)
plt.plot(ds, fitfunc(ds,*popt3), ls='-',lw=1.5,c='k')
plt.errorbar(densitys/1e18, phases[:,5], yerr=phases[:,6], ls="", marker='o',lw=1.5,c='k')



with io.open('cond_phase_vs_density.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()