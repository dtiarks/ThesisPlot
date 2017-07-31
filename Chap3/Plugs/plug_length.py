# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pandas as pd
import json
import io

def linear(x, a, b):
    return a*x + b
    
a0=c.physical_constants["Bohr radius"][0]
kB=c.k
m=87*c.physical_constants["atomic mass constant"][0]
lD1 = 780*10**-9

A2nd=(3.28)

profile1=np.loadtxt("img_004/profileavg.had")

#cut=np.loadtxt("img_004/Long_cut.had")
#ref_cut=np.loadtxt("img_006/cut.had")
ref_avg=np.loadtxt("img_002/profileavg.had")

bckg_list= np.empty(220)
bckg_list.fill(455.0)
prof1=profile1[:,1]-ref_avg[:,1]


frequencies = np.array([2.5,3.0,3.5,4.0,5.0,6.0,8.0,10.0,15.0])
lengths = np.array([78,97,118,138,179,220,304,386,597])

popt,pcov = opt.curve_fit(linear, frequencies, lengths)
print popt
print pcov
textstr = '$Slope=%.0f\mu m/MHz$\n$Offset=%.0f\mu m$'%(popt[0], popt[1])

plot_dict={}

plot1, = plt.plot(frequencies, lengths, 'bo', markersize = 10)
plot2, = plt.plot(np.arange(0.0,16.0,0.01), linear(np.arange(0.0,16.0, 0.01), popt[0], popt[1]), 'g', linewidth = 2.0)
plt.xlabel('frequency difference, $f_1 - f_2$ (MHz)', fontsize = 18)
plt.ylabel('FWHM ($\mu m$)', fontsize = 18)
plt.legend((plot1, plot2), ('data', 'fit'), fontsize = 20, loc = 'best', numpoints = 1)
#plt.grid(True, linewidth = 2)
#plt.title('Length vs. Frequency', fontsize = 20)
plt.tick_params(axis='both', which='major', labelsize=16, width = 2, size = 5)
plt.rcParams['axes.linewidth'] = 2
plt.ylim(0,650)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.text(10, 135, textstr, fontsize=18, verticalalignment='top', bbox=props)

h=pd.DataFrame(index=frequencies,data=lengths)
h2=pd.DataFrame(index=np.arange(0.0,16.0,0.01),data=linear(np.arange(0.0,16.0, 0.01), popt[0], popt[1]))
plot_dict['121']={
    'A':{'type':'plot','y':h2[0].to_json()} ,
    'B':{'type':'scatter','y':h[0].to_json(),'ylabel':r'FWHM ($\mu m$)','xlabel':r'Frequenzdifferenz $\Delta f_{AOD}$ (MHz)','label':'Daten','xlim':(0,16),'num':'a'}
                                   
}


h=pd.DataFrame(index=A2nd*profile1[:,0],data=prof1)
#h2=pd.DataFrame(index=np.arange(0.0,16.0,0.01),data=linear(np.arange(0.0,16.0, 0.01), popt[0], popt[1]))
plot_dict['122']={
    'A':{'type':'plot','y':h[0].to_json(),'xlabel':r'z-Koordinate ($\mu m$)','ylabel':r'OD (arb. einh.)','num':'b','xlim':(0,300)}
                                   
}

with io.open('pluglength.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))

#plt.tight_layout()
plt.show()
#plt.savefig('FreqVsLength_Poster.png', transparent = True)