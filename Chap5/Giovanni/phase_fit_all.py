import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
import re
from scipy.signal import butter, lfilter, firwin
from scipy.optimize import curve_fit

phaseLst=[]
phaseLstRef=[]
ampLst=[]
ampLstRef=[]

n=21

    
def fitFunc(t, A, phi,C):
    return A*np.sin(2*np.pi*20*t+phi)+C


for i in range(1,n+1):
   name="Phase_vs_fProbe.tsv_%03d"%i

   histo=np.loadtxt(("histo_%s")%name,skiprows=2000)
   histoRef=np.loadtxt(("histoRef_%s")%name,skiprows=2000)
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

freqFitP=np.loadtxt("phase-fit-control-off.txt")[:,0]
pfitOff=np.loadtxt("phase-fit-control-off.txt")[:,1]
pfitOn=np.loadtxt("phase-fit-control-on.txt")[:,1]

freqFitO=np.loadtxt("OD-fit-control-off.txt")[:,0]
OfitOff=np.loadtxt("OD-fit-control-off.txt")[:,1]
OfitOn=np.loadtxt("OD-fit-control-on.txt")[:,1]


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


from matplotlib import rc
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}
rc('font', **font)
fig=plt.figure(1, figsize=(8, 9), dpi=80,)
fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)

ax0 = fig.add_subplot(211)
ax0.clear()
ds=np.linspace(4,16,num=n)
dsf=np.linspace(4,16,num=len(freqFitO))
#ax0.set_title("$|r>=68 s_{1/2}$; $\Delta=10 MHz$$")
ax0.plot(ds,2*np.log(aRef[:,0]/a[:,0]),'ro',label="$\Omega_c=$???")
ax0.plot(ds,2*np.log(aRef[:,1]/a[:,1]),'bo',label="$\Omega_c=$0")
ax0.plot(dsf,OfitOff,"b")
ax0.plot(dsf,OfitOn,"r")

#ax0.plot(freq,aRef[:,0],'ro--',label="control on")
#ax0.plot(freq,aRef[:,1],'bo--',label="control off")
#ax0.set_xlabel("Probe detuning (MHz)")
ax0.tick_params(axis="x",which="both",labelbottom="off")
ax0.set_ylabel("OD")
#ax0.grid(True)
trans = ax0.get_xaxis_transform() # x in data untis, y in axes fraction
ax0.annotate('(a)', xy=(4.3,0.9 ), xycoords=trans)

#legds = ax0.legend(fancybox=True, framealpha=0.0,loc=1)

ax1 = fig.add_subplot(212)
ax1.clear()


dsf=np.linspace(4,16,num=len(freqFitP))
ax1.plot(ds,(p[:,0]-pRef[:,0]),'ro',label="control on")
ax1.plot(ds,(p[:,1]-pRef[:,1]),'bo',label="control off")
ax1.plot(dsf,pfitOff,"b")
ax1.plot(dsf,pfitOn,"r")

#ax1.plot(freq,pRef[:,0]/np.pi,'ro--',label="control on")
#ax1.plot(freq,pRef[:,1]/np.pi,'bo--',label="control off")
ax1.set_xlabel("Signal detuning (MHz)")
ax1.set_ylabel("Phase (rad)")
#dp = (popt0[1]-popt1[1])/np.pi
trans = ax1.get_xaxis_transform() # x in data untis, y in axes fraction
ax1.annotate('(b)', xy=(4.3,0.9 ), xycoords=trans)
#ax1.legend()
#ax1.grid(True)
y=p[:,0]-pRef[:,0]
y1=p[:,1]-pRef[:,1]
resadjust(ax1,yres=0.2)
plt.show()
fig.set_size_inches(7,6)
fig.savefig('phase_amplitude2.pdf')


