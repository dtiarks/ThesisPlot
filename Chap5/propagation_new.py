# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 17:28:05 2017

@author: daniel
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import io
import scipy.optimize as opt
from scipy import special
from scipy.optimize import curve_fit


def fitfunc(x, A):
    return A*x

def pol_komp(t,A,r_f,r_b, wd,d):
    Ip=1-special.erf((t-(wd+d))/r_f)*special.erf((t-d)/r_b)
    return A*Ip  
    


plot_dict={}

f=plt.figure()

propagation_dir="propagation_new"

histo = np.loadtxt(os.path.join(propagation_dir, "histo_delay_007.tsv"))

r_f_plus=0.2
r_b_plus=0.1
w_plus=3.
d_plus=11.

popt1=[0.001, r_f_plus, r_b_plus, w_plus, d_plus]
popt1, pcov1 = curve_fit(pol_komp, histo[:,0], histo[:,1], popt1)

popt2=[0.002, r_f_plus, r_b_plus, w_plus, d_plus]
popt2, pcov2 = curve_fit(pol_komp, histo[:,0], histo[:,2], popt2)

#h=pd.DataFrame(index=erf1[:,0],data=erf1[:,1])
#h2=pd.DataFrame(index=erf2[:,0],data=erf2[:,1])
#h3=pd.DataFrame(index=s_params[:,0],data=s_params[:,7])
#h3_err=pd.DataFrame(index=s_params[:,0],data=np.sqrt(s_params[:,7]))
#h4=pd.DataFrame(index=s_params[:,0],data=s_params[:,8])
#h4_err=pd.DataFrame(index=s_params[:,0],data=np.sqrt(s_params[:,8]))
#
#plot_dict['121']={
#    'A':{'type':'plot','y':h[0].to_json()},
#    'B':{'type':'plot','y':h2[0].to_json()},
#    'C':{'type':'errorbar','y':h3[0].to_json(),'yerr':h3_err[0].to_json(),'ylabel':u'Ereignisse','xlabel':u'Zeit ($\mu s$)','num':'a','xlim':(11.,14.5),'ylim':(0,80),'label':'Daten'},
#    'D':{'type':'errorbar','y':h4[0].to_json(),'yerr':h4_err[0].to_json()}
#}


txs=np.linspace(11.,14.5,200)


plt.subplot(121)

plt.plot(txs, pol_komp(txs,*popt1), ls='-',lw=1.5,c='r')
plt.plot(txs, pol_komp(txs,*popt2), ls='-',lw=1.5,c='b')
plt.scatter(histo[:,0], histo[:,1], marker='o',c='r')
plt.scatter(histo[:,0], histo[:,2], marker='o',c='b')
plt.xlim((11.,14.5))


I_p=pol_komp(txs,*popt1)
I_m=pol_komp(txs,*popt2)

S3=np.true_divide(I_p-I_m,I_p+I_m+1e-12)

theta=np.arctan(np.clip(1./(S3),-1e12,1e12))


I_p_data=histo[:,1]
I_m_data=histo[:,2]

S3_data=np.true_divide(I_p_data-I_m_data,I_p_data+I_m_data+1e-12)

theta_data=np.arctan(np.clip(1./(S3_data),-1e12,1e12))



#h=pd.DataFrame(index=polar[:-2,0],data=polar[:-2,1])
#h_err=pd.DataFrame(index=polar[:-2,0],data=polar[:-2,2])
#h2=pd.DataFrame(index=polar_fit[:,0],data=polar_fit[:,1])
#
#plot_dict['122']={
#    'A':{'type':'plot','y':h2[0].to_json()},
#    'B':{'type':'errorbar','y':h[0].to_json(),'yerr':h_err[0].to_json(),'ylabel':u'Polarwinkel (rad)','xlabel':u'Zeit ($\mu s$)','num':'b','ylim':(-1.5,1.5),'xlim':(11.,14.5),'label':'Daten'},
#}
#
plt.subplot(122)
plt.plot(txs, theta, ls='-',lw=1.5,c='r')
plt.scatter(histo[:,0], theta_data,c='r')
plt.xlim((11.,14.5))


#with io.open('propagation.json', 'w+') as f:
#  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()