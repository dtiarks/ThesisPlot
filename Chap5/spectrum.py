#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 14:51:43 2017

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
import re

phaseLst=[]
phaseLstRef=[]
ampLst=[]
ampLstRef=[]

n=21

    
def fitFunc(t, A, phi,C):
    return A*np.sin(2*np.pi*20*t+phi)+C


for i in range(1,n+1):
   name="Phase_vs_fProbe.tsv_%03d"%i

   histo=np.loadtxt(os.path.join("Giovanni","histo_%s"%name),skiprows=2000)
   histoRef=np.loadtxt(os.path.join("Giovanni","histoRef_%s"%name),skiprows=2000)
   histoInterval=(2530-2000,2580-2000)

   hT=histo[histoInterval[0]:histoInterval[1],0]
   hCrop1 = histo[histoInterval[0]:histoInterval[1],1]
   hCrop2 = histo[histoInterval[0]:histoInterval[1],2]
   
   hTRef=histoRef[histoInterval[0]:histoInterval[1],0]
   hCrop1Ref = histoRef[histoInterval[0]:histoInterval[1],1]
   hCrop2Ref = histoRef[histoInterval[0]:histoInterval[1],2]

   popt0=[0.01, 4,0.002]
   popt1=[0.01, 4,0.002]

   popt0, pcov0=curve_fit(fitFunc, hT, hCrop1, popt0)
   popt1, pcov0=curve_fit(fitFunc, hT, hCrop2, popt1)
   
   popt0Ref, pcov0Ref=curve_fit(fitFunc, hTRef, hCrop1Ref, popt0)
   popt1Ref, pcov0Ref=curve_fit(fitFunc, hTRef, hCrop2Ref, popt1)

   phaseLst.append([popt0[1],popt1[1]])
   ampLst.append([popt0[0],popt1[0]])
   phaseLstRef.append([popt0Ref[1],popt1Ref[1]])
   ampLstRef.append([popt0Ref[0],popt1Ref[0]])
   print("sample "+str(i)+" analyzed")

p=np.array(phaseLst)
a=np.array(ampLst)

pRef=np.array(phaseLstRef)
aRef=np.array(ampLstRef)

freq=np.linspace(-5,5,num=n)

freqFitP=np.loadtxt(os.path.join("Giovanni","phase-fit-control-off.txt"))[:,0]
pfitOff=np.loadtxt(os.path.join("Giovanni","phase-fit-control-off.txt"))[:,1]
pfitOn=np.loadtxt(os.path.join("Giovanni","phase-fit-control-on.txt"))[:,1]

freqFitO=np.loadtxt(os.path.join("Giovanni","OD-fit-control-off.txt"))[:,0]
OfitOff=np.loadtxt(os.path.join("Giovanni","OD-fit-control-off.txt"))[:,1]
OfitOn=np.loadtxt(os.path.join("Giovanni","OD-fit-control-on.txt"))[:,1]



fig=plt.figure(1, figsize=(8, 9), dpi=80,)
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)

ax0 = fig.add_subplot(211)
ax0.clear()
ds=np.linspace(4,16,num=n)
dsf=np.linspace(4,16,num=len(freqFitO))
ax0.plot(ds,2*np.log(aRef[:,0]/a[:,0]),'ro',label="$\Omega_c=$???")
ax0.plot(ds,2*np.log(aRef[:,1]/a[:,1]),'bo',label="$\Omega_c=$0")
ax0.plot(dsf,OfitOff,"b")
ax0.plot(dsf,OfitOn,"r")

ax0.tick_params(axis="x",which="both",labelbottom="off")
ax0.set_ylabel("OD")
trans = ax0.get_xaxis_transform() # x in data untis, y in axes fraction
ax0.annotate('(a)', xy=(4.3,0.9 ), xycoords=trans)


ax1 = fig.add_subplot(212)
ax1.clear()


dsf=np.linspace(4,16,num=len(freqFitP))
ax1.plot(ds,(p[:,0]-pRef[:,0]),'ro',label="control on")
ax1.plot(ds,(p[:,1]-pRef[:,1]),'bo',label="control off")
ax1.plot(dsf,pfitOff,"b")
ax1.plot(dsf,pfitOn,"r")

plt.show()

#plot_dict={}
#
#h=pd.DataFrame(index=ds,data=2*np.log(aRef[:,0]/a[:,0]))
#h2=pd.DataFrame(index=ds,data=2*np.log(aRef[:,0]/a[:,1]))
#h3=pd.DataFrame(index=dsf,data=OfitOff)
#h4=pd.DataFrame(index=dsf,data=OfitOn)
#plot_dict['121']={
#    'A':{'type':'scatter','y':h[0].to_json(),'num':'a','label':'Kontroll an','xlabel':u'Signalverstimmung $\Delta_s/2 \pi$ (MHz)','ylabel':u'OD','xlim':(4,16),'ylim':(0,3.2)},
#    'B':{'type':'scatter','y':h2[0].to_json(),'label':'Kontroll aus'},
#    'C':{'type':'plot','y':h3[0].to_json()},
#    'D':{'type':'plot','y':h4[0].to_json()},
#}
#
#dsf=np.linspace(4,16,num=len(freqFitP))
#
#h=pd.DataFrame(index=ds,data=(p[:,0]-pRef[:,0]),)
#h2=pd.DataFrame(index=ds,data=(p[:,1]-pRef[:,1]))
#h3=pd.DataFrame(index=dsf,data=pfitOff)
#h4=pd.DataFrame(index=dsf,data=pfitOn)
#plot_dict['122']={
#    'A':{'type':'scatter','y':h[0].to_json(),'num':'b','xlabel':u'Signalverstimmung $\Delta_s/2 \pi$ (MHz)','ylabel':u'Phase (rad)','xlim':(4,16)},
#    'B':{'type':'scatter','y':h2[0].to_json()},
#    'C':{'type':'plot','y':h3[0].to_json()},
#    'D':{'type':'plot','y':h4[0].to_json()},
#}
#
#with io.open('spectrum.json', 'w+') as f:
#  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))
