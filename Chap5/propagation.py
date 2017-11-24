# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 19:34:09 2017

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

propagation_dir="propagation"

s_params = np.loadtxt(os.path.join(propagation_dir, "S_parameters_Delay007.tsv"))
polar_fit = np.loadtxt(os.path.join(propagation_dir, "PolarAngleFit.txt"))
polar = np.loadtxt(os.path.join(propagation_dir, "PolarAngle.txt"),skiprows=233)
erf1 = np.loadtxt(os.path.join(propagation_dir, "Erf-fit1.txt"))
erf2 = np.loadtxt(os.path.join(propagation_dir, "Erf-fit2.txt"))




#h3_err=pd.DataFrame(index=densitys/1e18,data=phases[:,2])
#
#h4=pd.DataFrame(index=densitys/1e18,data=phases[:,3])
#h4_err=pd.DataFrame(index=densitys/1e18,data=phases[:,4])
#
#
#popt, pcov = opt.curve_fit(fitfunc, densitys/1e18, phases[:,1])
#popt2, pcov2 = opt.curve_fit(fitfunc, densitys/1e18, phases[:,3])
#popt3, pcov3 = opt.curve_fit(fitfunc, densitys/1e18, phases[:,5])
#
#ds=np.linspace(0,2.2,100)
#
h=pd.DataFrame(index=erf1[:,0],data=erf1[:,1])
h2=pd.DataFrame(index=erf2[:,0],data=erf2[:,1])
h3=pd.DataFrame(index=s_params[:,0],data=s_params[:,7])
h3_err=pd.DataFrame(index=s_params[:,0],data=np.sqrt(s_params[:,7]))
h4=pd.DataFrame(index=s_params[:,0],data=s_params[:,8])
h4_err=pd.DataFrame(index=s_params[:,0],data=np.sqrt(s_params[:,8]))

plot_dict['121']={
    'A':{'type':'plot','y':h[0].to_json()},
    'B':{'type':'plot','y':h2[0].to_json()},
    'C':{'type':'errorbar','y':h3[0].to_json(),'yerr':h3_err[0].to_json(),'ylabel':u'Ereignisse','xlabel':u'Zeit ($\mu s$)','num':'a','xlim':(11.,14.5),'ylim':(0,80),'label':'Daten'},
    'D':{'type':'errorbar','y':h4[0].to_json(),'yerr':h4_err[0].to_json()}
}

plt.subplot(121)

plt.plot(erf1[:,0], erf1[:,1], ls='-',lw=1.5,c='r')
plt.plot(erf2[:,0], erf2[:,1], ls='-',lw=1.5,c='b')
plt.scatter(s_params[:,0], s_params[:,7], marker='o',c='r')
plt.scatter(s_params[:,0], s_params[:,8], marker='o',c='b')
plt.xlim((11.,14.5))



h=pd.DataFrame(index=polar[:-2,0],data=polar[:-2,1])
h_err=pd.DataFrame(index=polar[:-2,0],data=polar[:-2,2])
h2=pd.DataFrame(index=polar_fit[:,0],data=polar_fit[:,1])

plot_dict['122']={
    'A':{'type':'plot','y':h2[0].to_json()},
    'B':{'type':'errorbar','y':h[0].to_json(),'yerr':h_err[0].to_json(),'ylabel':u'Polarwinkel (rad)','xlabel':u'Zeit ($\mu s$)','num':'b','ylim':(-1.5,1.5),'xlim':(11.,14.5),'label':'Daten'},
}

plt.subplot(122)
plt.plot(polar_fit[:,0], polar_fit[:,1], ls='-',lw=1.5,c='r')
plt.errorbar(polar[:,0], polar[:,1], yerr=polar[:,2], ls="", marker='o',lw=1.5,c='k')
plt.xlim((11.,14.5))


with io.open('propagation.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()