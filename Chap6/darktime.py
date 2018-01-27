# -*- coding: utf-8 -*-
"""
Created on Sat Jan 27 11:57:11 2018

@author: daniel
"""

import locale
# Set to German locale to get comma decimal separater
locale.setlocale(locale.LC_NUMERIC, 'deu_deu')
import Tomography as tom
#import Tomography_old_detection as tom_o
import EIT
import quPy as qp
#from scipy.optimize import curve_fit
import numpy as np
import matplotlib as mpl
import pylab as plt
import pandas as pd
import json
import io
from scipy.optimize import curve_fit
#mpl.use('pgf')

LC = '7' # LoadCycles
data_dir = './darktime_data/'#'W:/data/2017/01/260117/'

def simple_unwrap(phases):
    unwrapped=[phases[0]]
    for p in phases[1:]:
        d = (p - unwrapped[-1]) % (2*np.pi) # distance between phases, restricted to 2 pi
        d = d - (np.sign(d-np.pi)+1)/2. * 2*np.pi # shift upper half downwards
        print d
        unwrapped.append(unwrapped[-1] + d)
    return unwrapped

filelist = [data_dir + "lifetime_001/photons%d" % i for i in xrange(1, 2100)]
datafields = [qp.MicrosecondsSinceTrigger]
categorizers = [qp.POL, qp.ROI, qp.N_EXP]
filters={'TomographyBase':['2']}#{'ROI':[0]} # we only want to look at R and L
t_min=0
t_max=1000
from_fileheader=['TomographyBases', 'QCompCh2Delay', 'QCompCh4Delay', 'CycleID', 'LoadCycles']
profile = [datafields, categorizers, filters, t_min, t_max, from_fileheader]
ROIs = None

T = tom.Tomography(filelist, profile, ROIs)

DELAYS = list(T.get_categories()['QCompCh2Delay'])
DELAYS.sort(lambda x,y: cmp(float(x),float(y)))

DARKTIMES = [float(d) - float(DELAYS[0]) + 0.05 for d in DELAYS] # first scan step was 200 ns dark time

ROIs = [[10,11.05]]+[[11.05+d, 12.05+d] for d in DARKTIMES]

T = tom.Tomography(filelist, profile, ROIs)
             
# get reference clickrates (without atoms)
c_ref, c_ref_err = T.getCounts(['pol','ROI','LoadCycles'], True)

ryd_ref = c_ref[('L',0,'0')]
ryd_ref_err = c_ref_err[('L',0,'0')]
gs_ref = c_ref[('R',0,'0')]
gs_ref_err = c_ref_err[('R',0,'0')]
                       
dop, doperr, pol, polerr, azi, azierr = T.getVectorCounts(['pol','ROI','LoadCycles', 'QCompCh2Delay'])
vis, viserr = T.getEquatorialDOPCounts(['pol','ROI','LoadCycles', 'QCompCh2Delay'])
dop0, doperr0, pol0, polerr0, azi0, azierr0 = [el[(0,'0')] for el in T.getVectorCounts(['pol','ROI','LoadCycles'])]

print "Input: \nPolarwinkel %.2f +- %.2f\nAzimut %.2f +- %.2f\nPhotonenzahl %.1f" % (pol0, polerr0, azi0, azierr0, (ryd_ref+gs_ref)*10)

RYD_EFFS = []
RYD_EFF_ERRS = []
GS_EFFS = []
GS_EFF_ERRS = []
AZIs = []
AZI_ERRS = []

c, c_err = T.getCounts(['pol','ROI', 'LoadCycles', 'QCompCh2Delay'], True)

i=1
for d in DELAYS:
    RYD_EFFS.append(c[('L',i,LC,d)]/ryd_ref)
    RYD_EFF_ERRS.append( RYD_EFFS[-1] * np.sqrt((c_err[('L',i,LC,d)]/c[('L',i,LC,d)])**2 + (ryd_ref_err/ryd_ref)**2) )
    GS_EFFS.append(c[('R',i,LC,d)]/gs_ref)
    GS_EFF_ERRS.append( GS_EFFS[-1] * np.sqrt((c_err[('R',i,LC,d)]/c[('R',i,LC,d)])**2 + (gs_ref_err/gs_ref)**2) )
    AZIs.append(azi[(i,LC,d)])
    AZI_ERRS.append(azierr[(i,LC,d)])
    i+=1
    
# try to unwrap the angles and do a linear fit
AZIs=np.array(AZIs)
DARKTIMES = np.array(DARKTIMES)
AZIs=np.concatenate((AZIs[:5], AZIs[5:10], AZIs[10:11]+2*np.pi ) )

f = lambda t, phi0, dpdt: phi0 + 2*np.pi*dpdt*t
g = lambda t, alpha, deltarm: np.angle(1 + alpha*np.exp(-1j*deltarm*t))
h = lambda t, phi0, dpdt, alpha, deltarm: f(t, phi0, dpdt) + g(t, alpha, deltarm)
eta = lambda t, eta0, tau, alpha, deltarm: np.exp(-t/tau) * eta0 * (1 + alpha**2 + 2 * alpha * np.cos(deltarm*t)) / (1+alpha)**2
popt, pcov = curve_fit(f , DARKTIMES, AZIs, p0 = [-2, 2*np.pi*0.15])
AZI2 = AZIs - f(DARKTIMES, *popt)
#popt2, pcov2 = curve_fit(g, DARKTIMES, AZI2, [0.9, 1.4])
#popt, pcov = curve_fit(h , DARKTIMES, AZIs, p0 = [-2, 2*np.pi*0.15, 0.55,  1.36])
popt0, pcov0 = curve_fit(eta, DARKTIMES[1:], RYD_EFFS[1:], p0 = [0.2, 2, 0.58,  1.3])
x_fit = np.linspace(0, 8, 500)

#fit_desc = r"$\varphi (t_D) = \SI{%s}{} + 2\pi t_D \cdot \SI{%s}{\mega\hertz}$" % tuple([un2str(*pair) for pair in zip(popt, np.sqrt(np.diag(pcov)))])
fit_desc = ""


plt.figure(0)
ax = plt.subplot(211)
plt.ylabel('Effizienz')
plt.errorbar(DARKTIMES, GS_EFFS, yerr=GS_EFF_ERRS, marker='o', color='orange', label='R')#, ls='none')
plt.errorbar(DARKTIMES, RYD_EFFS, yerr=RYD_EFF_ERRS, marker='o', color = '#85c1fe', label='L')#, ls='none')
#plt.plot(x_fit, eta(x_fit, *popt0), color='#85c1fe')
plt.legend(title = 'Polarisation')
#plt.margins(0.1,0.1)
plt.ylim(0,0.2)
plt.xlim(0,7.5)
ax.set_xticklabels([])

#np.savetxt("gs_effs.tsv", np.column_stack((DARKTIMES, GS_EFFS, GS_EFF_ERRS)))
#np.savetxt("ryd_effs.tsv", np.column_stack((DARKTIMES, RYD_EFFS, RYD_EFF_ERRS)))

plt.subplot(212)
plt.xlabel(r'Dunkelzeit $t_D$ ($\si{\micro\second}$)')
plt.ylabel('Azimut (rad)')
plt.plot(x_fit, f(x_fit, *popt), 'k-', color = 'red', label=fit_desc)
#plt.plot(x_fit, f(x_fit, *popt) + g(x_fit, *popt2), 'g-')
plt.errorbar(DARKTIMES, AZIs, yerr=AZI_ERRS, marker='o', color='red', ls='none')
plt.xlim(0,7.5)
plt.ylim(-3.5,4.8)
plt.legend(fancybox = True, loc='upper left')

#np.savetxt("azi.tsv", np.column_stack((DARKTIMES, AZIs, AZI_ERRS)))

plt.tight_layout()
plt.show()
#savefig("dark_time_scan")