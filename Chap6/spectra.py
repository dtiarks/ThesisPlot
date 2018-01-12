# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 12:36:31 2017

@author: mpqqubec
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
mpl.use('pgf')

RYD_BLUE = '#85c1fe'  # the color for the rydberg data

BINNING=10

# First Rydberg EIT
LC='5'
data_dir = './spectra_data/'

files = [data_dir+"eit_sweep_001/photons%d" % i for i in xrange(1,65)]
datafields = [qp.MillisecondsSinceStart]
categorizers = [qp.POL, qp.N_EXP, qp.ROI]
from_fileheader = ["LoadCycles", "CycleID", "TomographyBases"]
profile = [datafields, categorizers, {}, 0, 1000, from_fileheader]

S = EIT.EITSpectrum(files, profile, [[11,15],[20,100]])

c, cerr = S.getCounts(level=['LoadCycles','pol','ROI'], per_exp = True)

print "Input Photon Number: ca. %.1f |L>, %.1f |R>" % (c[('0','L',0)]*10, c[('0','R',0)]*10)

trans, trans_err = S.getTransmission(BINNING, [0,1000], ['LoadCycles', 'pol', 'ROI'], 
                                     signal = ('5','R',0), 
                                     reference = ('0','R',0), 
                                     fit_reference = False,
                                     use_fit_errors = False,
                                     fit_order = 1,
                                     background = None)#('5','R',1))


# initial parameters for the fit
S.OD = 30
S.Oc = 10 # all in MHz
S.Dc = 0
S.gp = 0.4
S.f0 = 0
S.G = 6.07
S.loss = np.inf # this means no loss

y = trans.values
yErr = trans_err.values
x = (trans.index.values-500)*60/1000 # in MHz
  
popt, perr, f = S.fitSpectrum('EIT+OD+f0', x, y, yErr)
center = 500 + S.f0 / 60. * 1000
x = (trans.index.values-center)*60./1000.  # recalibrate x-axis with fitted f0 parameter
popt, perr, f = S.fitSpectrum('EIT+OD+f0', x, y, yErr)
x_fit = np.linspace(-30,30,300)
y_fit = f(x_fit, *popt)

#desc_string = r'''$\Delta_c/2\pi = \SI{%s}{\mega\hertz}$
#$\gamma_d/2\pi = \SI{%s}{\micro\second\tothe{-1}}$
#$\Omega_c/2\pi = \SI{%s}{\mega\hertz}$
#OD$_{\mathrm{max}} = \SI{%s}{}$'''  % tuple([un2str(*pair) for pair in zip(popt[:4], perr[:4])])

desc_string = ''


T = tom.Tomography(files, profile, [[11,15]])
dop, dop_err, pol, pol_err, azi, azi_err = T.getVectorHisto(BINNING, [0,1000], 'abs', ['LoadCycles', 'ROI'])
#azi = azi*(-1) # azimut has the wrong sign because we exchanged R and L in the setup, BUT: we want to plot a phaseshift, not the azimut rotation
offset = azi[('0',0)].mean()

print "Phase offset: %.2f +- %.2f" % (offset, np.std(azi[('0',0)])/np.sqrt(len(azi[('0',0)])))

# unwrap the phase data
phase = lambda X: S.XnRe(X,S.Dc,S.gp,S.Oc,S.f0,28,S.loss,S.G) # use a OD of 28 for unwrapping, gives better results
azi_ref = np.transpose([phase(x) for n in xrange(-3,4)]) # create a reference array of the same size as the one below, containing the fit values 
azi_choices = np.transpose([azi[(LC,0)].values - offset + n*np.pi*2 for n in xrange(-3, 4)]) # several arrays containing 2*pi up- and downshifted versions of the datapoints
idx = np.argmin(np.abs(azi_choices-azi_ref), axis=1) # determine the indices of the elements with the smallest deviation from the theory curve
azi_unwrapped = [azi_choices[i][idx[i]] for i in xrange(len(azi_choices))] # retrieve those elements and pack them into a new array

phase = lambda X: S.XnRe(X,S.Dc,S.gp,S.Oc,S.f0,S.OD,S.loss,S.G) # for plotting use the fitted OD

#f = plt.figure(0)
#ax = plt.subplot(211)
#plt.plot(x_fit, y_fit, color='#85c1fe')
#plt.errorbar(x,y, yerr=yErr, marker='o', ls='none', color='#85c1fe')
#plt.ylabel("Transmission")
#plt.ylim(0,1)
#plt.xlim(-30,30)
#ax.set_xticklabels([])
#plt.text(3,0.4, desc_string)
#plt.subplot(212)
#plt.plot(x_fit, phase(x_fit) , color='#85c1fe')
#plt.errorbar(x, azi_unwrapped, yerr=azi_err[(LC,0)].values, marker='o', ls='none', color='#85c1fe')
#plt.xlabel(r"Verstimmung $ \Delta_s / 2\pi$ $(MHz)$")
#plt.ylabel(r"$\Delta \phi$ (rad)")
#plt.xlim(-30,30)


plot_dict={}

h=pd.DataFrame(index=x_fit,data=y_fit)
h2=pd.DataFrame(index=x,data=y)
h2_err=pd.DataFrame(index=x,data=yErr)

plot_dict['221']={
    'A':{'type':'plot','y':h[0].to_json(),'ylabel':r'Transmission','xlim':(-30,30),'ylim':(0,1),'num':'a'} ,
    'B':{'type':'errorbar','y':h2[0].to_json(),'yerr':h2_err[0].to_json()}
                                   
}

h=pd.DataFrame(index=x_fit,data=phase(x_fit))
h2=pd.DataFrame(index=x,data=azi_unwrapped)
h2_err=pd.DataFrame(index=x,data=azi_err[(LC,0)].values)

plot_dict['223']={
    'A':{'type':'plot','y':h[0].to_json(),'xlabel':r"Verstimmung $ \Delta_s / 2\pi$ (MHz)",'ylabel':r"$\Delta \phi$ (rad)",'xlim':(-30,30),'ylim':(-9.7,9.7),'num':'c'} ,
    'B':{'type':'errorbar','y':h2[0].to_json(),'yerr':h2_err[0].to_json()}
                                   
}






# Now Groundstate EIT
LC='7'
data_dir = './spectra_data/'

files = [data_dir+"eit_002/photons%d" % i for i in xrange(1,80)]

S = EIT.EITSpectrum(files, profile, [[11,15],[20,100]])

c, cerr = S.getCounts(level=['LoadCycles','pol','ROI'], per_exp = True)

print "Input Photon Number: ca. %.1f |L>, %.1f |R>" % (c[('0','L',0)]*10, c[('0','R',0)]*10)

trans, trans_err = S.getTransmission(BINNING, [0,1000], ['LoadCycles', 'pol', 'ROI'], 
                                     signal = ('7','L',0), 
                                     reference = ('0','L',0), 
                                     fit_reference = False,
                                     use_fit_errors = False,
                                     fit_order = 1,
                                     background = None)#('5','R',1))

S.structure='lambda'
# initial parameters for the fit
S.OD = 6
S.Oc = 5 # all in MHz
S.Dc = 0
S.gp = 0.4
S.f0 = 0
S.G = 6.07
S.loss = np.inf # this means no loss

y = trans.values
yErr = trans_err.values
x = (trans.index.values-500)*-20/1000 # in MHz, using the -1. AOM diffraction order
  
popt, perr, f = S.fitSpectrum('EIT+OD+f0', x, y, yErr)
center = 500 + S.f0 / -20. * 1000
x = (trans.index.values-center)*-20./1000.  # recalibrate x-axis with fitted f0 parameter
popt, perr, f = S.fitSpectrum('EIT+OD+f0', x, y, yErr)
x_fit = np.linspace(-10,10,300)
y_fit = f(x_fit, *popt)

#desc_string = r'''$\Delta_c/2\pi = \SI{%.3g \pm %.1g}{\mega\hertz}$
#$\gamma_d/2\pi = \SI{%.3g \pm %.1g}{\micro\second\tothe{-1}}$
#$\Omega_c/2\pi = \SI{%.3g \pm %.1g}{\mega\hertz}$
#OD$_{\mathrm{max}} = \SI{%.3g \pm %.1g}{}$'''  % tuple(np.ravel( zip(popt[:4], perr[:4])))
#desc_string = r'''$\Delta_c/2\pi = \SI{%s}{\mega\hertz}$
#$\gamma_d/2\pi = \SI{%s}{\micro\second\tothe{-1}}$
#$\Omega_c/2\pi = \SI{%s}{\mega\hertz}$
#OD$_{\mathrm{max}} = \SI{%s}{}$'''  % tuple([un2str(*pair) for pair in zip(popt[:4], perr[:4])])

T = tom.Tomography(files, profile, [[11,15]])
dop, dop_err, pol, pol_err, azi, azi_err = T.getVectorHisto(BINNING, [0,1000], 'abs', ['LoadCycles', 'ROI'])
azi = azi*(-1) # azimut has the wrong sign because we exchanged R and L in the setup
offset = azi[('0',0)].mean()

print "Phase offset: %.2f +- %.2f" % (offset, np.std(azi[('0',0)])/np.sqrt(len(azi[('0',0)])))

# unwrap the phase data
phase = lambda X: S.XnRe(X,S.Dc,S.gp,S.Oc,S.f0,S.OD,S.loss,S.G) 
azi_ref = np.transpose([phase(x) for n in xrange(-3,4)]) # create a reference array of the same size as the one below, containing the fit values 
azi_choices = np.transpose([azi[(LC,0)].values - offset + n*np.pi*2 for n in xrange(-3, 4)]) # several arrays containing 2*pi up- and downshifted versions of the datapoints
idx = np.argmin(np.abs(azi_choices-azi_ref), axis=1) # determine the indices of the elements with the smallest deviation from the theory curve
azi_unwrapped = [azi_choices[i][idx[i]] for i in xrange(len(azi_choices))] # retrieve those elements and pack them into a new array
azi_unwrapped_errs = list(azi_err[(LC,0)].values)
x_unwrapped = list(x)

i=0
while i < len(azi_unwrapped_errs): # remove all points with too large errorbars
    if azi_unwrapped_errs[i] > np.pi/2:
        del azi_unwrapped_errs[i]
        del azi_unwrapped[i]
        del x_unwrapped[i]
    else:
        i+=1

#f = plt.figure(1)
#ax = plt.subplot(211)
#plt.plot(x_fit, y_fit, color='orange')
#plt.errorbar(x,y, yerr=yErr, marker='o', ls='none', color='orange')
#plt.ylabel("Transmission")
#plt.ylim(0,1)
#plt.xlim(-10,10)
#ax.set_xticklabels([])
#plt.text(1,0.3, desc_string)
#plt.subplot(212)
#plt.plot(x_fit, phase(x_fit) , color='orange')
#plt.errorbar(x_unwrapped, azi_unwrapped, yerr=azi_unwrapped_errs, marker='o', color='orange', ls='none')
#plt.xlabel(r"Verstimmung $ \Delta_s / 2\pi (MHz)$")
#plt.ylabel(r"$\Delta \phi$ (rad)")
#plt.xlim(-10,10)


h=pd.DataFrame(index=x_fit,data=y_fit)
h2=pd.DataFrame(index=x,data=y)
h2_err=pd.DataFrame(index=x,data=yErr)

plot_dict['222']={
    'A':{'type':'plot','y':h[0].to_json(),'ylabel':r'Transmission','xlim':(-10,10),'ylim':(0,1),'num':'b'} ,
    'B':{'type':'errorbar','y':h2[0].to_json(),'yerr':h2_err[0].to_json()}
                                   
}

h=pd.DataFrame(index=x_fit,data=phase(x_fit))
h2=pd.DataFrame(index=x_unwrapped,data=azi_unwrapped)
h2_err=pd.DataFrame(index=x_unwrapped,data=azi_unwrapped_errs)

plot_dict['224']={
    'A':{'type':'plot','y':h[0].to_json(),'xlabel':r"Verstimmung $ \Delta_s / 2\pi$ (MHz)",'ylabel':r"$\Delta \phi$ (rad)",'xlim':(-10,10),'num':'d'} ,
    'B':{'type':'errorbar','y':h2[0].to_json(),'yerr':h2_err[0].to_json()}
                                   
}



with io.open('memory_spectra.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))

#plt.show()
