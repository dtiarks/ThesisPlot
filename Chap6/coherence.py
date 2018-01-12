# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 17:55:53 2018

@author: daniel
"""

import numpy as np
import matplotlib as mpl
import pylab as plt
import pandas as pd
import json
import io
import os

data_dir = "./coherence_data"

azi = np.loadtxt(os.path.join(data_dir, "memory_Azimuth_120117_header.txt"))
t_azi = np.loadtxt(os.path.join(data_dir, "azi_time.txt"))
vis = np.loadtxt(os.path.join(data_dir, "memory_Visibility_120117_header.txt"))
counts_lr = np.loadtxt(os.path.join(data_dir, "memory_CountsLR_120117_header.txt"))
counts_l_fit = np.loadtxt(os.path.join(data_dir, "Counts-L-fit.txt"))
counts_r_fit = np.loadtxt(os.path.join(data_dir, "Counts-R-fit.txt"))

t = np.loadtxt(os.path.join(data_dir, "time_corrected.txt"))

watt_factor = 2.54591

counts_lr[:,1:] = counts_lr[:,1:]*watt_factor
counts_l_fit [:,1:] = counts_l_fit [:,1:]*watt_factor
counts_r_fit [:,1:] = counts_r_fit [:,1:]*watt_factor


f = plt.figure(0)

ax = plt.subplot(311)

plt.errorbar(t,counts_lr[:,1],yerr=counts_lr[:,2], c='b', label='R', ls='', marker='o')
plt.plot(counts_r_fit [:,0],counts_r_fit [:,1])

plt.errorbar(t,counts_lr[:,3],yerr=counts_lr[:,4], c='r', label='L', ls='', marker='o')
plt.plot(counts_l_fit [:,0],counts_l_fit [:,1])
plt.xlim(-0.1,0.9)

plot_dict={}


h_L=pd.DataFrame(index=t,data=counts_lr[:,1])
h_L_err=pd.DataFrame(index=t,data=counts_lr[:,2])
h_R=pd.DataFrame(index=t,data=counts_lr[:,3])
h_R_err=pd.DataFrame(index=t,data=counts_lr[:,4])

h_L_fit=pd.DataFrame(index=counts_r_fit [:,0],data=counts_l_fit [:,1])
h_R_fit=pd.DataFrame(index=counts_r_fit [:,0],data=counts_r_fit [:,1])

plot_dict['311']={
    'A':{'type':'errorbar','y':h_L[0].to_json(),'yerr':h_L_err[0].to_json(),'ylabel':r'Leistung (fW)','num':'a','label':'L','xlim':(-0.1,0.9),'ylim':(0,40)} ,
    'B':{'type':'errorbar','y':h_R[0].to_json(),'yerr':h_R_err[0].to_json(),'label':'R'},
    'C':{'type':'plot','y':h_L_fit[0].to_json()},
    'D':{'type':'plot','y':h_R_fit[0].to_json()}
                                   
}


ax = plt.subplot(312)

plt.errorbar(t_azi,azi[:,1],yerr=azi[:,2], c='b', label='R', ls='', marker='o')
plt.xlim(-0.1,0.9)
plt.ylim(0,2*np.pi)

h=pd.DataFrame(index=t_azi[1:7],data=azi[1:7,1])
h_err=pd.DataFrame(index=t_azi[1:7],data=azi[1:7,2])

plot_dict['312']={
    'A':{'type':'errorbar','y':h[0].to_json(),'yerr':h_err[0].to_json(),'ylabel':r'Azimuth (rad)','num':'b','xlim':(-0.1,0.9),'ylim':(0,2*np.pi)}                                    
}

ax = plt.subplot(313)

plt.errorbar(t_azi,vis[:,1],yerr=vis[:,2], c='b', label='R', ls='', marker='o')
plt.xlim(-0.1,0.9)
plt.ylim(0,1)

h=pd.DataFrame(index=t_azi[1:7],data=vis[1:7,1])
h_err=pd.DataFrame(index=t_azi[1:7],data=vis[1:7,2])

plot_dict['313']={
    'A':{'type':'errorbar','y':h[0].to_json(),'yerr':h_err[0].to_json(),'ylabel':r'Visibility','xlabel':r'Zeit ($\mu s$)','num':'c','xlim':(-0.1,0.9),'ylim':(0,1)}                                    
}

with io.open('memory_coherence.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))

plt.show()

