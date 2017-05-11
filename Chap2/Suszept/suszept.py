# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 11:26:50 2017

@author: daniel
"""

import Tomography as tom
import quPy as qp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import io


def susceptibility(Delta_s, Delta_c, gamma_21, Omega_c, ladder=True):
    delta = (Delta_s + (-1 + 2*int(ladder)) * Delta_c)
    return 1j*(gamma_21 - 2j * delta)/(np.abs(Omega_c)**2 + (1 - 2j * Delta_s)*(gamma_21 - 2j * delta))

# the parameters, all in units of Gamma_3
Delta_c = 0
gamma_21 = 0
Omega_c = 1.5

# calculate the curve
detuning = np.linspace(-3, 3, 400)
chi = susceptibility(detuning, Delta_c, gamma_21, Omega_c)
chi_2L = susceptibility(detuning, Delta_c, gamma_21, 0)

plot_dict={}

# plot it
f=plt.figure()
plt.subplot(221)
plt.plot(detuning, np.imag(chi), ls='-',lw=1.5,color='b')
plt.plot(detuning, np.imag(chi_2L), ls='-',lw=1.5,color='g')
plt.axvline(0, ls='--',lw=1.0)
#plt.xlabel(u'$\Delta_s / \Gamma_3$')
plt.ylabel(u'Im($\chi / \chi_0$)')
plt.margins(0,0.1)
#plt.tight_layout()

h=pd.DataFrame(index=detuning,data=np.imag(chi))
h2=pd.DataFrame(index=detuning,data=np.imag(chi_2L))
plot_dict['221']={
    'A':{'type':'plot','y':h[0].to_json(),'ylabel':u'Im($\chi / \chi_0$)','margin':(0,0.0),'num':'a'},                
    'B':{'type':'plot','y':h2[0].to_json(),'margin':(0,0.0)},
    'C':{'type':'axv','y':0}
}

plt.subplot(223)
plt.plot(detuning, np.real(chi), ls='-',lw=1.5)
plt.plot(detuning, np.real(chi_2L), ls='-',lw=1.5)
#plt.xlabel(u'$\Delta_s / \Gamma_3$')
plt.ylabel(u'Re($\chi / \chi_0$)')
plt.yticks(np.arange(-0.4, 0.6, 0.2))
plt.margins(0,0.1)
plt.axvline(0, ls='--',lw=1.0)
#plt.tight_layout()

h=pd.DataFrame(index=detuning,data=np.real(chi))
h2=pd.DataFrame(index=detuning,data=np.real(chi_2L))
plot_dict['223']={
    'A':{'type':'plot','y':h[0].to_json(),'ylabel':u'Re($\chi / \chi_0$)','xlabel':u'$\Delta_s / \Gamma_3$','margin':(0,0.0),'num':'c'},                
    'B':{'type':'plot','y':h2[0].to_json(),'margin':(0,0.0)},
    'C':{'type':'axv','y':0}
}

# detuned eit
# the parameters, all in units of Gamma_3
Delta_c = 1.5
gamma_21 = 0
Omega_c = 1.5

# calculate the curve
detuning = np.linspace(-3, 3, 400)
chi = susceptibility(detuning, Delta_c, gamma_21, Omega_c)
chi_2L = susceptibility(detuning, Delta_c, gamma_21, 0)


plt.subplot(222)
plt.plot(detuning, np.imag(chi), ls='-',lw=1.5)
plt.plot(detuning, np.imag(chi_2L), ls='-',lw=1.5)
plt.xlabel(u'$\Delta_s / \Gamma_3$')
#plt.ylabel(u'Im($\chi / \chi_0$)')
plt.margins(0,0.1)
plt.axvline(-1.5, ls='--',lw=1.0)
#plt.tight_layout()

h=pd.DataFrame(index=detuning,data=np.imag(chi))
h2=pd.DataFrame(index=detuning,data=np.imag(chi_2L))
plot_dict['222']={
    'A':{'type':'plot','y':h[0].to_json(),'margin':(0,0.0),'num':'b'},                
    'B':{'type':'plot','y':h2[0].to_json(),'margin':(0,0.0)},
    'C':{'type':'axv','y':-1.5}
}

plt.subplot(224)
plt.plot(detuning, np.real(chi), ls='-',lw=1.5)
plt.plot(detuning, np.real(chi_2L), ls='-',lw=1.5)
plt.xlabel(u'$\Delta_s / \Gamma_3$')
#plt.ylabel(u'Re($\chi / \chi_0$)')
plt.yticks(np.arange(-0.4, 0.6, 0.2))
plt.margins(0,0.1)
plt.axvline(-1.5, ls='--',lw=1.0)
#plt.tight_layout()

h=pd.DataFrame(index=detuning,data=np.real(chi))
h2=pd.DataFrame(index=detuning,data=np.real(chi_2L))
plot_dict['224']={
    'A':{'type':'plot','y':h[0].to_json(),'xlabel':u'$\Delta_s / \Gamma_3$','margin':(0,0.0),'num':'d'},                
    'B':{'type':'plot','y':h2[0].to_json(),'margin':(0,0.0)},
    'C':{'type':'axv','y':-1.5}
}

with io.open('suszept.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()