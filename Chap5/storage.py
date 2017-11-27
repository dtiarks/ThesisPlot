# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 19:17:43 2017

@author: daniel
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import io

def decay(t,A,tau,f,V):
  return A*np.exp(-t/tau)*(1 - V*np.sin(2*np.pi*f*t)**2)
#*np.exp(-t/tau2)  
def decay2(t,A,tau):
  return A*np.exp(-2*t**2/tau**2)
  
def trilobite_decay(t,alpha,omega,tau,A0):
    return A0*np.exp(-t**1/(tau**1))*(np.cos(alpha)**4+np.sin(alpha)**4+2*np.cos(alpha)**2*np.sin(alpha)**2*np.cos(omega*t))


tauTh=11.7

store_dir="store_retr"
  
data1=np.loadtxt(os.path.join(store_dir,"store_001.tsv"))

data1[:,0]=data1[:,0]+0.15

#popt1=[0.2,5,0.25,0.4]
popt1=[0.5,2*np.pi*0.2,5.,0.2]
popt1, pcov1 = curve_fit(trilobite_decay, data1[:,0],data1[:,1], popt1)

print 1e3*popt1[1]/(2.*np.pi)
print 2.*np.pi/popt1[1]

#fitRange = np.arange(0,15,0.01)
#plt.plot(data1[:,0],(data1[:,1]),"o--",label=r'Data (store_001)')
#plt.plot(fitRange,decay(fitRange, popt1[0], popt1[1], popt1[2], popt1[3]))

txs=np.linspace(0,16,100)

plot_dict={}

f=plt.figure()


h=pd.DataFrame(index=data1[:,0],data=data1[:,1])
h_err=pd.DataFrame(index=data1[:,0],data=data1[:,2])
h2=pd.DataFrame(index=txs,data=trilobite_decay(txs,*popt1))


plot_dict['111']={
    'A':{'type':'errorbar','y':h[0].to_json(),'yerr':h_err[0].to_json(),'ylabel':u'Effizienz','xlabel':u'Speicherzeit ($\mu s$)','xlim':(0.,12.),'ylim':(0,0.25),'label':'Daten'},
    'B':{'type':'plot','y':h[0].to_json()}
}

plt.subplot(111)
plt.errorbar(data1[:,0], data1[:,1], yerr=data1[:,2], ls="", marker='o',lw=1.5,c='r')
plt.plot(txs, trilobite_decay(txs,*popt1))


#h=pd.DataFrame(index=densitys/1e18,data=phases[:,5])
#h_err=pd.DataFrame(index=densitys/1e18,data=phases[:,6])
#h2=pd.DataFrame(index=ds,data=fitfunc(ds,*popt3))

#plot_dict['122']={
#    'A':{'type':'plot','y':h2[0].to_json()},
#    'B':{'type':'errorbar','y':h[0].to_json(),'yerr':h_err[0].to_json(),'xlabel':u'Dichte ($10^{12}\ cm^{-3}$)','num':'b','xlim':(0.,2.5),'ylim':(0,4.0),'label':'Daten'},
#}

#plt.subplot(122)
#plt.plot(ds, fitfunc(ds,*popt3), ls='-',lw=1.5,c='k')
#plt.errorbar(densitys/1e18, phases[:,5], yerr=phases[:,6], ls="", marker='o',lw=1.5,c='k')



#with io.open('storage_retrieval.json', 'w+') as f:
#  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()