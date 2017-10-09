# -*- coding: utf-8 -*-
"""
Created on Mon Oct 09 13:29:17 2017

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
    
popt, pcov = opt.curve_fit(fitfunc, Temp, freq)
print popt

plot_dict={}

f=plt.figure()


ring_down = np.loadtxt("cavity_data.dat", skiprows=2)
data_len = len(ring_down[:,1])
data_idx=np.arange(data_len)
data_sample = np.sort(np.random.choice(data_idx,5000, replace=False))

t_list=np.array(ring_down[65500:,1])

popt_rd, pcov_rd = opt.curve_fit(decay, t_list, ring_down[65500:,2],p0=[13.1, 6.9, 0.1, -0.04])

ts=np.linspace(t_list[0],t_list[-1],1000)
fitted_rd=decay(ts, popt_rd[0], popt_rd[1], popt_rd[2], popt_rd[3])

rd_time=np.array(ring_down[data_sample,1])
rd_signal=np.array(ring_down[data_sample,2])

h=pd.DataFrame(index=rd_time,data=rd_signal)
h2=pd.DataFrame(index=ts,data=fitted_rd)
plot_dict['121']={
    'A':{'type':'plot','y':h[0].to_json(),'ylabel':u'Photodiodensignal (mV)','xlabel':u'Zeit ($\mu s$)','num':'a','xlim':(0.,50.),'ylim':(-5,35),'label':'Daten'},
    'B':{'type':'plot','y':h2[0].to_json(),'label':'Kurvenanpassung'}
}

plt.subplot(121)
plt.plot(ring_down[data_sample,1] , ring_down[data_sample,2], ls='-',lw=1.5,c='r')
plt.plot(ts , fitted_rd, ls='-',lw=1.5,c='b')


fitted=fitfunc(Ts, popt[0], popt[1], popt[2])
h=pd.DataFrame(index=Temp,data=freq)
h2=pd.DataFrame(index=Ts,data=fitted)

plot_dict['122']={
    'A':{'type':'scatter','y':h[0].to_json(),'xlabel':u'Temperatursollwert ($^\circ\mathrm{C}$)','num':'b','ylabel':u'Frequenzdifferenz (MHz)'},
    'B':{'type':'plot','y':h2[0].to_json()}
}

plt.subplot(122)
plt.scatter(Temp, freq ,c='r')
plt.plot(Ts, fitted, ls='-',lw=1.5,c='b')



with io.open('cavity_characterization.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()

