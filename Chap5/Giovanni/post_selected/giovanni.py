import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
import re
from scipy.optimize import curve_fit
from matplotlib import rc

#curve=np.loadtxt("curve.tsv")

histoLst=[]
histoPostLst=[]
histoRefLst=[]

xLst=np.arange(-8,-12.5,-1)
#print xLst

name="ref001"
nameRef="ref001"

##for n in ["histo_%03d.tsv"%i for i in range(2,5)]:
##  a=np.loadtxt(n)
##  #print a
##  histoLst.append(np.asarray(a))
##
##
##for n in ["histoPost_%03d.tsv"%i for i in range(2,5)]:
##  histoPostLst.append(np.loadtxt(n))
##
##for n in ["histoRef_%03d.tsv"%i for i in range(2,5)]:
##  histoRefLst.append(np.loadtxt(n))

histo=np.loadtxt("histo%s.tsv"%name)
histoPost=np.loadtxt("histoPost%s.tsv"%name)
histoRef=np.loadtxt("histoRef%s.tsv"%name)
meanRef=np.loadtxt("meanRef%s.tsv"%name)
numEvents=np.loadtxt("numEvents%s.tsv"%name)
ifile=open("histo%s.tsv"%name, "r")
s=ifile.read()
RoiTarget=tuple(float(x)*1e6 for x in re.search("'RoiTarget', \((.*?),(.*?)\)", s).groups())
RoiRetrieval=tuple(float(x)*1e6 for x in re.search("'RoiRetrieval', \((.*?),(.*?)\)", s).groups())
ifile.close()


histoRef2=np.loadtxt("histo%s.tsv"%nameRef)
histoPostRef=np.loadtxt("histoPost%s.tsv"%nameRef)
histoRefRef=np.loadtxt("histoRef%s.tsv"%nameRef)
meanRefRef=np.loadtxt("meanRef%s.tsv"%nameRef)
ifileRef=open("histo%s.tsv"%nameRef, "r")
sRef=ifileRef.read()
RoiTargetRef=tuple(float(x)*1e6 for x in re.search("'RoiTarget', \((.*?),(.*?)\)", sRef).groups())
RoiRetrievalRef=tuple(float(x)*1e6 for x in re.search("'RoiRetrieval', \((.*?),(.*?)\)", sRef).groups())
ifileRef.close()

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

pLst=[]
pErrLst=[]

pLstRef=[]
pErrLstRef=[]

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 13}
rc('font', **font)

#for step in range(0,5):
for step in [0]:
  fig=plt.figure(0)
  fig2=plt.figure(1)
  fig3=plt.figure(2)
  fig.clear()

  ax1=fig2.add_subplot(111)
#  ax12=fig2.add_subplot(212)
  ax2=fig.add_subplot(311)
  ax3=fig.add_subplot(312)
  ax4=fig.add_subplot(313)
  fig2.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=0.15)
  fig3.subplots_adjust(left=None, bottom=0.2, right=None, top=None, wspace=None, hspace=0.15)
  axR=fig3.add_subplot(111)
   
    #  ax1.set_title("Step %d:" % step+ " Evaluated " +datetime.now().strftime("%c\n")
    #  + "Incoming: $N_g$=%.2f " % meanRef[step,0]
    #  + "$N_t$=%.2f " % meanRef[step,4]
#  )
    
  
  ax1.plot(histo[:,step*6]-delay, histo[:,step*6+2], color='r',label='Gate on')
  axR.plot(histoRef[:,step*6]-delay, histoRef[:,step*6+2], color='g', ls="-")
  axR.set_xlim([0,12.55-delay])
  axR.set_ylim([0,0.11])
#  axR.set_ylabel("transmitted photons\n in 50ns",labelpad=-0.05)
#  ax12.plot(histo[:,step*6]-delay, histo[:,step*6+4], color='b',label='Gate off')
#  ax12.plot(histoRef[:,step*6]-delay, histoRef[:,step*6+4], color='b', ls="--")
  resadjust(ax1,yres=0.02,xres=0.3)
  resadjust(axR,yres=0.04,xres=0.3)
#  resadjust(ax12,yres=0.02,xres=0.3)
  ax1.set_xlim([0,12.55-delay])
  ax1.set_ylim([0,0.054])
  ax1.set_ylabel("transmitted photons\n in 50ns",labelpad=-0.05)
  ax1.set_xlabel("Time ($\mu s$)",labelpad=0)
#  ax1.legend()
  trans = ax1.get_xaxis_transform()
#  ax1.annotate('(b)', xy=(0.1,.85 ), xycoords=trans)
#  trans = ax12.get_xaxis_transform()
#  ax12.annotate('(b)', xy=(0.1,.85 ), xycoords=trans)
  axR.tick_params(axis="x",which="both",labelbottom="off")

#  ax12.set_xlim([0,12.55-delay])
#  ax12.set_ylim([0,0.08])
#  ax12.set_ylabel("transmitted photons\n in 50ns",labelpad=-0.05)
  
#  ax12.legend()

  b1 = int(histoInterval[0]/(binningHisto*10**6))
  b2 = int(histoInterval[1]/(binningHisto*10**6))
  
  hT1=histoPost[b1:b2,step*9]
  hT2=histo[b1:b2,step*6]
  hCrop1 = histoPost[b1:b2,step*9+2]/histoPost[b1:b2,step*9+4]*3.7
  #hCrop1 = histoPost[b1:b2,step*9+2]
  
  hCrop2 = histo[b1:b2,step*6+4]
  
  
  popt0=fitParameters["initPost"]
  popt1=fitParameters["initRef"]

  popt0, pcov0=curve_fit(fitFunc, hT1, hCrop1, popt0)
  popt1, pcov1=curve_fit(fitFunc, hT2, hCrop2, popt1)
  
  #popt0=fitParameters["initPost"]
  #popt1=fitParameters["initRef"]
  #print popt0
  #print popt1

  hCropRef1 = histoRef[b1:b2,step*6+2]
  hCropRef2 = histoRef[b1:b2,step*6+4]

  poptRef0=fitParameters["initRef"]
  poptRef1=fitParameters["initRef"]

  poptRef0, pcovRef0=curve_fit(fitFunc, hT2, hCropRef1, poptRef0)
  poptRef1, pcovRef1=curve_fit(fitFunc, hT2, hCropRef2, poptRef1)
  
  fT1=np.linspace(hT1[0],hT1[-1],1e3)
  fT2=np.linspace(hT2[0],hT2[-1],1e3)
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
#  ax2.grid(True)
#  ax2.legend()
  resadjust(ax4,yres=0.01)

  ax3.plot(histo[:,step*6]-delay, histo[:,step*6+4], color='b',ls='', marker='o',label='No gate')
  ax3.plot(fT2-delay,fitFunc(fT2,*popt1),'b-')
  ax3.set_xlim([fitParameters["xscale"][0]-delay,fitParameters["xscale"][1]-delay])
  ax3.set_ylim([0,0.04])
  ax3.plot([1.5878,1.5878],[0,0.08],"k--")
  ax3.text(0.02,0.8,'(b)', transform=ax3.transAxes)
#  ax3.grid(True)
#  ax3.set_xlabel("Time ($\mu s$)")
  ax3.tick_params(axis="x",which="both",labelbottom="off")
#  ax3.tick_params(axis="x",which="both",labelbottom="off")
#  ax3.legend()
  resadjust(ax3,yres=0.01)
  
  ax2.plot(histoRef[:,step*6]-delay, histoRef[:,step*6+2], color='g',ls='', marker='o',label='w/o Atoms')
  ax2.plot(fT2-delay,fitFunc(fT2,*poptRef0),'g-')
  ax2.set_xlim([fitParameters["xscale"][0]-delay,fitParameters["xscale"][1]-delay])
  ax2.set_ylim([0,0.08])
  ax2.plot([1.5878,1.5878],[0,0.08],"k--")
  ax2.text(0.02,0.8,'(a)', transform=ax2.transAxes)
#  ax4.grid(True)
  ax4.set_xlabel("Time ($\mu s$)")
#  ax4.legend()
    
#  ax2.set_ylabel("transmitted\n photons\n in 50ns",labelpad=-0.05)
  ax3.set_ylabel("transmitted photons in 50ns",labelpad=-0.00, fontsize=14)    
#  ax4.set_ylabel("transmitted photons\n in 50ns",labelpad=-0.05)  
  resadjust(ax2,yres=0.02)
  
  p1=popt0[1]-poptRef0[1]
  p1Err=np.sqrt(pcov0[1,1]+pcovRef0[1,1])
  p2=popt1[1]-poptRef0[1]
  p2Err=np.sqrt(pcov1[1,1]+pcovRef0[1,1])
  dp = (p2-p1)
  dpErr = np.sqrt(p1Err**2+p2Err**2)
  #dpErr = np.sqrt(p1Err**2+p2Err**2)
  OD1=2*np.log(poptRef0[0]/popt0[0])
  OD2=2*np.log(poptRef0[0]/popt1[0])
  ODer1=np.sqrt((2./poptRef0[0])**2*pcovRef0[0][0]+(2./popt0[0])**2*pcov0[0][0])
  ODer2=np.sqrt((2./poptRef0[0])**2*pcovRef0[0][0]+(2./popt1[0])**2*pcov1[0][0])
  
  
  pLst.append(dp)
  pErrLst.append(dpErr)

  print "Postselected phase shift in pi"
  print "phase on "+str(p1)+" +- "+str(p1Err)
  print "phase off "+str(p2)+" +- "+str(p2Err)
  print "phase difference "+ str(dp) +" +- "+str(dpErr)
  
  print "Od on "+str(OD1)+" +- "+str(ODer1)
  print "OD off "+str(OD2)+" +- "+str(ODer2)
  
  print "w/o atoms"
  print poptRef1
  print "Gate off (full ensemble)"
  print popt1
  print "Postselected (gate on)"
  print popt0


plt.show()
fig.set_size_inches(10,6)
fig2.set_size_inches(6.3,2.5)
fig3.set_size_inches(6.3,2.5)
#fig2.savefig("transmitted0.pdf")
#fig.savefig("GateOnOff.pdf")
#fig3.savefig("reference.pdf")