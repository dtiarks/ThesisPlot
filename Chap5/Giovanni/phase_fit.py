import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
import re
from scipy.signal import butter, lfilter, firwin
from scipy.optimize import curve_fit
from matplotlib import rc

name="Phase_vs_fProbe.tsv_011"

photonsPerClick=3.7

histo=np.loadtxt(("histo_%s")%name)
histoInterval=(2530,2580)
#histoInterval=(2200,2880)

histoRef=np.loadtxt(("histoRef_%s")%name)

eit=np.loadtxt(("eit_%s")%name)


    
def fitFunc(t, A, phi,C):
    return A*np.sin(2*np.pi*20*t+phi)+C


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}
rc('font', **font)
#fig=plt.figure(1, figsize=(8, 12), dpi=80,)
fig=plt.figure(1,figsize=(7, 4.5))
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)


ax0 = fig.add_subplot(111)
ax0.clear()

hT=histo[histoInterval[0]:histoInterval[1],0]
hCrop1 = histo[histoInterval[0]:histoInterval[1],1]
hCrop2 = histo[histoInterval[0]:histoInterval[1],2]
hCropref = histoRef[histoInterval[0]:histoInterval[1],2]

popt0=[0.005, 1.35,0.002]
popt1=[0.005, 1.35,0.002]

popt0, pcov0=curve_fit(fitFunc, hT, hCrop1, popt0)
popt1, pcov1=curve_fit(fitFunc, hT, hCrop2, popt1)
poptR, pcovR=curve_fit(fitFunc, hT, hCropref, popt1)
fitT=np.linspace(hT[0],hT[-1],1e3)

#ax0.plot(np.asarray(hT)-12.5,hCrop1*photonsPerClick,'ro',label="$\Omega_c=$???")
#ax0.plot(np.asarray(fitT)-12.5,fitFunc(fitT,*popt0)*photonsPerClick,'r-')
# ax0.plot(np.asarray(hT)-12.5,hCrop2*photonsPerClick,'bo',label="$\Omega_c=$0")
# ax0.plot(np.asarray(fitT)-12.5,fitFunc(fitT,*popt1)*photonsPerClick,'b-')
# ax0.plot(np.asarray(hT)-12.5,hCropref*photonsPerClick,'go',label="$\Omega_c=$0")
# ax0.plot(np.asarray(fitT)-12.5,fitFunc(fitT,*poptR)*photonsPerClick,'g-')

# ax0.set_xlabel("Time ($\mu s$)")
# ax0.set_ylabel("Number of transmitted photons in 50ns")
# print "Phase of : without atoms, without control, with control "+str([poptR[1],popt1[1],popt0[1]])
# print "errors: without atoms, without control, with control "+str([np.sqrt(pcovR[1][1]),np.sqrt(pcov1[1][1]),np.sqrt(pcov0[1][1])])
# dpOff = (popt1[1]-poptR[1])
# dpOn = (popt0[1]-poptR[1])
# print dpOff,dpOn, dpOff-dpOn
# ErrOff=np.sqrt(pcov1[1][1]+pcovR[1][1])
# print ErrOff
# #leg = ax0.legend(fancybox=True, framealpha=0.0,loc=1)

# #ax0.set_title(("Amplitude modulation: $|r>=68 s_{1/2}; \phi=%.3f \pi$")%dp)
# #ax0.legend()
# #ax0.grid(True)



#ax1 = fig.add_subplot(12)
#ax1.clear()

#x=eit[100:120,0]
#freq=np.linspace(-10,10,num=len(x))

#y1=eit[100:120,1]
#y2=eit[100:120,2]

#ax1.plot(freq,y1,'r-',label="control on")
#ax1.plot(freq,y2,'b-',label="control off")
#ax1.set_title("Transmission")
#ax1.set_ylabel("Amplitude")
#ax1.set_xlabel("Probe detuning (MHz)")
#ax1.grid(True)

fig.set_size_inches(7,4.5)
plt.show()
#fig.set_size_inches(7,12)
# plt.savefig('timetrace_d0_1.pdf')


datOn=histo[2300:2750,1]
t=histo[2300:2750,0]
datOff=histo[2300:2750,2]
datRef=histoRef[2300:2750,1]
poptOn,pcovOn=curve_fit(fitFunc,t,datOn,[0,0,0])
poptOff,pcovOff=curve_fit(fitFunc,t,datOff,[0,0,0])
poptRef,pcovRef=curve_fit(fitFunc,t,datRef,[0,0,0])




