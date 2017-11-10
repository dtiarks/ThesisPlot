#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 11:01:51 2017

@author: daniel
"""

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
import re


def fitfunc(T, T0, A, B):
    return A*(T - T0)**2 + B

def decay(t, t0, tau, A, B):
    return A*np.exp(-(t-t0)/tau)+B

def fitFunc(t, A, phi,C):
  return A*np.cos(2*np.pi*20*(t-11.25)+phi)+C

def resadjust(ax, xres=None, yres=None):
    """
    Send in an axis and I fix the resolution as desired.
    """

    if xres:
        start, stop = ax.get_xlim()
        ticks = np.arange(start, stop + xres, xres)
        ax.set_xticks(ticks)
    if yres:
        start, stop = ax.get_ylim()
        ticks = np.arange(start, stop + yres, yres)
        ax.set_yticks(ticks)


fitParameters={
  "binningHisto":0.05e-7,
  "histoInterval":(11.25,11.75),
  "refScale":1/0.0002,
  "initRef":[25, 2,25],
  "initPost":[25, 2,25],
  "xscale":[11.05,11.7],
  "delay":9.94
}




binningHisto=fitParameters["binningHisto"]
histoInterval=fitParameters["histoInterval"]
delay=fitParameters["delay"]
A=fitParameters["refScale"]

histoLst=[]
histoPostLst=[]
histoRefLst=[]

xLst=np.arange(-8,-12.5,-1)

name="ref001"
nameRef="ref001"


histo=np.loadtxt(os.path.join("Giovanni","post_selected","histo%s.tsv"%name))
histoPost=np.loadtxt(os.path.join("Giovanni","post_selected","histoPost%s.tsv"%name))
histoRef=np.loadtxt(os.path.join("Giovanni","post_selected","histoRef%s.tsv"%name))
meanRef=np.loadtxt(os.path.join("Giovanni","post_selected","meanRef%s.tsv"%name))
numEvents=np.loadtxt(os.path.join("Giovanni","post_selected","numEvents%s.tsv"%name))
ifile=open(os.path.join("Giovanni","post_selected","histo%s.tsv"%name), "r")
s=ifile.read()
RoiTarget=tuple(float(x)*1e6 for x in re.search("'RoiTarget', \((.*?),(.*?)\)", s).groups())
RoiRetrieval=tuple(float(x)*1e6 for x in re.search("'RoiRetrieval', \((.*?),(.*?)\)", s).groups())
ifile.close()

#popt, pcov = opt.curve_fit(fitfunc, Temp, freq)
#print popt

step=0

plot_dict={}

fig=plt.figure()

ax2=fig.add_subplot(311)
ax3=fig.add_subplot(312)
ax4=fig.add_subplot(313)

b1 = int(histoInterval[0]/(binningHisto*10**6))
b2 = int(histoInterval[1]/(binningHisto*10**6))

hT1=histoPost[b1:b2,step*9]
hT2=histo[b1:b2,step*6]
hCrop1 = histoPost[b1:b2,step*9+2]/histoPost[b1:b2,step*9+4]*3.7
#hCrop1 = histoPost[b1:b2,step*9+2]

hCrop2 = histo[b1:b2,step*6+4]


popt0=fitParameters["initPost"]
popt1=fitParameters["initRef"]

popt0, pcov0=opt.curve_fit(fitFunc, hT1, hCrop1, popt0)
popt1, pcov1=opt.curve_fit(fitFunc, hT2, hCrop2, popt1)

#popt0=fitParameters["initPost"]
#popt1=fitParameters["initRef"]
#print popt0
#print popt1

hCropRef1 = histoRef[b1:b2,step*6+2]
hCropRef2 = histoRef[b1:b2,step*6+4]

poptRef0=fitParameters["initRef"]
poptRef1=fitParameters["initRef"]

poptRef0, pcovRef0=opt.curve_fit(fitFunc, hT2, hCropRef1, poptRef0)
poptRef1, pcovRef1=opt.curve_fit(fitFunc, hT2, hCropRef2, poptRef1)

fT1=np.linspace(hT1[0],hT1[-1],200)
fT2=np.linspace(hT2[0],hT2[-1],200)
ax4.errorbar(histoPost[:,step*9]-delay, (histoPost[:,step*9+2]/histoPost[:,step*9+4])*3.7, 0*np.sqrt(histoPost[:,step*9+2])/histoPost[:,step*9+4] *3.7 , color='r', ls='', marker='o',label='Gate on')
  #ax2.errorbar(histoPost[:,step*9], np.asarray(histoPost[:,step*9+2])/4500*3.7, np.sqrt(np.asarray(histoPost[:,step*9+2]))/numEvents[step,3] , color='r', ls='', marker='o',label='Postselected (on)')
##  ax2.errorbar(histoPost[:,step*9], histoPost[:,step*9+2], np.sqrt(histoPost[:,step*9+2]), color='r', ls='', marker='o',label='Postselected (on)')
ax4.plot(fT1-delay,fitFunc(fT1,*popt0),'r-')
ax4.set_xlim([fitParameters["xscale"][0]-delay,fitParameters["xscale"][1]-delay])
ax4.set_ylim([0,0.04])
ax2.tick_params(axis="x",which="both",labelbottom="off")
trans = ax4.get_xaxis_transform()
ax4.text(0.02,0.8,'(c)', transform=ax4.transAxes)
ax4.plot([1.5878,1.5878],[0,0.08],"k--")

resadjust(ax4,yres=0.01)

ax3.plot(histo[:,step*6]-delay, histo[:,step*6+4], color='b',ls='', marker='o',label='No gate')
ax3.plot(fT2-delay,fitFunc(fT2,*popt1),'b-')
ax3.set_xlim([fitParameters["xscale"][0]-delay,fitParameters["xscale"][1]-delay])
ax3.set_ylim([0,0.04])
ax3.plot([1.5878,1.5878],[0,0.08],"k--")
ax3.text(0.02,0.8,'(b)', transform=ax3.transAxes)

ax3.tick_params(axis="x",which="both",labelbottom="off")

resadjust(ax3,yres=0.01)

ax2.plot(histoRef[:,step*6]-delay, histoRef[:,step*6+2], color='g',ls='', marker='o',label='w/o Atoms')
ax2.plot(fT2-delay,fitFunc(fT2,*poptRef0),'g-')
ax2.set_xlim([fitParameters["xscale"][0]-delay,fitParameters["xscale"][1]-delay])
ax2.set_ylim([0,0.08])
ax2.plot([1.5878,1.5878],[0,0.08],"k--")
ax2.text(0.02,0.8,'(a)', transform=ax2.transAxes)

ax4.set_xlabel("Time ($\mu s$)")

ax3.set_ylabel("transmitted photons in 50ns",labelpad=-0.00, fontsize=14)

resadjust(ax2,yres=0.02)



h=pd.DataFrame(index=histo[2100:2500,step*6]-delay,data=histo[2100:2500,step*6+4])
h2=pd.DataFrame(index=fT2-delay,data=fitFunc(fT2,*popt1))
plot_dict['121']={
    'A':{'type':'scatter','y':h[0].to_json(),'num':'a','xlabel':u'Zeit ($\mu s$)',
         'xlim':(fitParameters["xscale"][0]-delay,fitParameters["xscale"][1]-delay),'ylim':(0,0.05),'label':'Daten'},
    'B':{'type':'plot','y':h2[0].to_json(),'label':'Kurvenanpassung','ylabel':u'Intensit\"at (Photonen/50 ns)'},
    'C':{'type':'axv','y':1.5878}
}


h=pd.DataFrame(index=histoPost[2100:2500,step*9]-delay,data=(histoPost[2100:2500,step*9+2]/histoPost[2100:2500,step*9+4])*3.7)
h2=pd.DataFrame(index=fT1-delay,data=fitFunc(fT1,*popt0))

plot_dict['122']={
    'A':{'type':'scatter','y':h[0].to_json(),'xlabel':u'Zeit ($\mu s$)',
         'xlim':(fitParameters["xscale"][0]-delay,fitParameters["xscale"][1]-delay),'num':'b','ylim':(0,0.05)},
    'B':{'type':'plot','y':h2[0].to_json()},
    'C':{'type':'axv','y':1.5878}
}

plt.show()

with io.open('sideband_postselected_phaseshift.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))




