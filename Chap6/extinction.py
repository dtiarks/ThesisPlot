# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 10:53:04 2017

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

watt_factor = 2.54591

BINNING = 0.05
LC = '6' # LoadCycles
data_dir = './extinction_data/'

filelist = [data_dir + "gs_store_001/photons%d" % i for i in xrange(1, 264)]
datafields = [qp.MicrosecondsSinceTrigger]
categorizers = [qp.POL, qp.SPCM, qp.ROI, qp.N_EXP]
t_min=0
t_max=1000
from_fileheader=['TomographyBases', 'CycleID', 'ExecExtParam0', 'LoadCycles']
profile = [datafields, categorizers, {}, t_min, t_max, from_fileheader]
ROIs = [[10.4,11.1],[15.2,16.4],[20,100]]
r0, r1, r2 = np.ravel(np.diff(ROIs, axis=1)) # ROI lengths, for count rates

# determine dark counts from all files
T = tom.Tomography(filelist, profile, ROIs)
c, c_err = T.getCounts(['LoadCycles','ROI','pol'],True)
dark_L = 0 #c[(LC,2,'L')]/r2 * 1e6 * 1e-3
dark_L_err = 0 # c_err[(LC,2,'L')]/r2 * 1e6 * 1e-3
dark_R = 0 # c[(LC,2,'R')]/r2 * 1e6 * 1e-3
dark_R_err = 0 #c_err[(LC,2,'R')]/r2 * 1e6 * 1e-3
print T.getCounts(['ROI','spcm'],True)[0][2]/r2 *1e6
print "Dark count rates: %.4f +- %.4f kHz L and %.4f +- %.4f kHz R" % (dark_L, dark_L_err, dark_R, dark_R_err)

profile[2]={'ExecExtParam0':'0.3'} # filter out lowest photon number
T = tom.Tomography(filelist, profile, ROIs)
h_ret, h_ret_err = T.getHisto(BINNING,[14.5,17.], 'rel', ['LoadCycles','pol'], True)
h_ret = h_ret[LC]/BINNING*1e6*1e-3
h_ret_err = h_ret_err[LC]/BINNING*1e6*1e-3
h_in, h_in_err =  T.getHisto(BINNING,[10.,11.5], 'rel', ['LoadCycles','pol'], True)
h_in = h_in['0']/BINNING*1e6*1e-3
h_in_err = h_in_err['0']/BINNING*1e6*1e-3
                   
# calculate extinction
c, c_err = T.getCounts(['LoadCycles','ROI','pol'],True)
cR = c[(LC,1,'R')]/r1 * 1e6 * 1e-3 # divide by ROI width in µs and multiply by 1/(µ*k) to get mean count rate in kHz
cRe = c_err[(LC,1,'R')]/r1 * 1e6 * 1e-3
cL = c[(LC,1,'L')]/r1 * 1e6 * 1e-3
cLe = c_err[(LC,1,'L')]/r1 * 1e6 * 1e-3
cR0 = c[('0',0,'R')]/r0 * 1e6 * 1e-3
cRe0 = c_err[('0',0,'R')]/r0 * 1e6 * 1e-3
cL0 = c[('0',0,'L')]/r0 * 1e6 * 1e-3
cLe0 = c_err[('0',0,'L')]/r0 * 1e6 * 1e-3
             
eta_rw = 100 * cR*r1 / (cR0*r0) # must multiply by ROI width again to get correct efficiency
eta_rw_err = eta_rw * np.sqrt((cRe/cR)**2 + (cRe0/cR0)**2)

N_in = (c[('0',0,'R')] + c[('0',0,'L')])*10
                 
print "Retrieval: %.2f +- %.2f kHz R %.2f +- %.2f kHz L\nInput: %.2f +- %.2f kHz R %.2f +- %.2f kHz L" % (cR, cRe, cL, cLe, cR0, cRe0, cL0, cLe0)

ext = (cL- dark_L)/(cR - dark_R)
ext_err = np.sqrt( (cLe/(cR-dark_R))**2 + (dark_L_err/(cR-dark_R))**2 + (ext*cRe/(cR-dark_R))**2 + (ext*dark_R_err/(cR - dark_R))**2 )
ext0 = (cL0 - dark_L)/(cR0- dark_R)
ext0_err = np.sqrt( (cLe0/(cR0-dark_R))**2 + (dark_L_err/(cR0-dark_R))**2 + (ext0*cRe0/(cR0-dark_R))**2 + (ext0*dark_R_err/(cR0 - dark_R))**2 )

bgc=(.8,.8,.8,1) # background color for the axes

f = plt.figure(0)
ax = plt.subplot(221, axisbg=bgc)
ax.axvspan(*ROIs[0], color='white')
plt.errorbar(h_in[('L')].index.values, h_in[('L')], yerr=h_in_err[('L')].values, marker = 'o', color='#85c1fe', label="L")
plt.errorbar(h_in[('R')].index.values, h_in[('R')], yerr=h_in_err[('R')].values, marker = 'o', color='orange', label="R")
plt.ylabel(u"Zählrate (kHz)")
ax.locator_params(tight=True, nbins=6)
plt.legend(title='Polarisation', fancybox=True, loc='upper left')
plt.margins(0,.1)
plt.ylim(ymin=0)
#plt.text(10.83,20, r'$N_{\mathrm{in}} = \SI{%.1f}{}$' % N_in + '\n' + 
#                  r'$\epsilon_{R\mathrm{, in}} = \SI{%s}{\percent}$' % un2str(ext0*100, ext0_err*100), ha='center' )

plot_dict={}

h_L=pd.DataFrame(index=h_in[('L')].index.values-11.0,data=watt_factor*h_in[('L')].values)
h_L_err=pd.DataFrame(index=h_in[('L')].index.values-11.0,data=watt_factor*h_in_err[('L')].values)
h_R=pd.DataFrame(index=h_in[('R')].index.values-11.0,data=watt_factor*h_in[('R')].values)
h_R_err=pd.DataFrame(index=h_in[('R')].index.values-11.0,data=watt_factor*h_in_err[('R')].values)

plot_dict['221']={
    'A':{'type':'errorbar','y':h_L[0].to_json(),'yerr':h_L_err[0].to_json(),'ylabel':r'Leistung (fW)','num':'a','label':'L','xlim':(-0.6,0.2001),'ylim':(-watt_factor*16,watt_factor*160)} ,
    'B':{'type':'errorbar','y':h_R[0].to_json(),'yerr':h_R_err[0].to_json(),'label':'R'}
                                   
}

ax = plt.subplot(222, axisbg=bgc)
ax.axvspan(*ROIs[1], color='white')
plt.errorbar(h_ret[('L')].index.values, h_ret[('L')], yerr=h_ret_err[('L')].values, marker='o', color='#85c1fe', label='L')
plt.errorbar(h_ret[('R')].index.values, h_ret[('R')], yerr=h_ret_err[('R')].values, marker='o', color='orange', label='R')
ax.locator_params(tight=True, nbins=6)
plt.margins(0,.1)
plt.ylim(ymin=0)
#plt.text(16,8,  r'$\epsilon_R = \SI{%s}{\percent}$' % (un2str(ext*100, ext_err*100)) + '\n' +
#                r'$\eta_{\mathrm{rw}} = \SI{%s}{\percent}$' % un2str(eta_rw, eta_rw_err), ha='center' )

h_L=pd.DataFrame(index=h_ret[('L')].index.values-11.0,data=watt_factor*h_ret[('L')].values)
h_L_err=pd.DataFrame(index=h_ret[('L')].index.values-11.0,data=watt_factor*h_ret_err[('L')].values)
h_R=pd.DataFrame(index=h_ret[('R')].index.values-11.0,data=watt_factor*h_ret[('R')].values)
h_R_err=pd.DataFrame(index=h_ret[('R')].index.values-11.0,data=watt_factor*h_ret_err[('R')].values)

plot_dict['222']={
    'A':{'type':'errorbar','y':h_L[0].to_json(),'yerr':h_L_err[0].to_json(),'num':'b','label':'L','xlim':(4.0,5.5),'ylim':(-watt_factor*1.4,watt_factor*14)} ,
    'B':{'type':'errorbar','y':h_R[0].to_json(),'yerr':h_R_err[0].to_json(),'label':'R'}
                                   
}

# LOWER PLOTS

BINNING = 0.03
LC = '6' # LoadCycles
#data_dir = 'W:/data/2017/01/260117/'
#data_dir = '../../../DataAnalysis/260117/'


filelist = [data_dir + "rydberg_store_001/photons%d" % i for i in xrange(1, 116)]
datafields = [qp.MicrosecondsSinceTrigger]
categorizers = [qp.POL, qp.SPCM, qp.ROI, qp.N_EXP]
t_min=0
t_max=1000
from_fileheader=['TomographyBases', 'CycleID', 'ExecExtParam0', 'LoadCycles']
profile = [datafields, categorizers, {}, t_min, t_max, from_fileheader]
ROIs = [[10.4,11.1],[15.2,15.9],[20,100]]
r0, r1, r2 = np.ravel(np.diff(ROIs, axis=1)) # ROI lengths, for count rates

# determine dark counts from all files
T = tom.Tomography(filelist, profile, ROIs)
c, c_err = T.getCounts(['LoadCycles','ROI','pol'],True)
dark_L = 0 #c[(LC,2,'L')]/r2 * 1e6 * 1e-3
dark_L_err = 0# c_err[(LC,2,'L')]/r2 * 1e6 * 1e-3
dark_R = 0 #c[(LC,2,'R')]/r2 * 1e6 * 1e-3
dark_R_err = 0 #c_err[(LC,2,'R')]/r2 * 1e6 * 1e-3
print T.getCounts(['ROI','spcm'],True)[0][2]/r2 *1e6
print "Dark count rates: %.4f +- %.4f kHz L and %.4f +- %.4f kHz R" % (dark_L, dark_L_err, dark_R, dark_R_err)

profile[2]={'ExecExtParam0':'0.3'} # filter out lowest photon number
T = tom.Tomography(filelist, profile, ROIs)
h_ret, h_ret_err = T.getHisto(BINNING,[14.5,17.], 'rel', ['LoadCycles','pol'], True)
h_ret = h_ret[LC]/BINNING*1e6*1e-3
h_ret_err = h_ret_err[LC]/BINNING*1e6*1e-3
h_in, h_in_err =  T.getHisto(BINNING,[10.,11.4], 'rel', ['LoadCycles','pol'], True)
h_in = h_in['0']/BINNING*1e6*1e-3
h_in_err = h_in_err['0']/BINNING*1e6*1e-3
                   
# calculate extinction
c, c_err = T.getCounts(['LoadCycles','ROI','pol'],True)
cR = c[(LC,1,'R')]/r1 * 1e6 * 1e-3 # divide by ROI width in µs and multiply by 1/(µ*k) to get mean count rate in kHz
cRe = c_err[(LC,1,'R')]/r1 * 1e6 * 1e-3
cL = c[(LC,1,'L')]/r1 * 1e6 * 1e-3
cLe = c_err[(LC,1,'L')]/r1 * 1e6 * 1e-3
cR0 = c[('0',0,'R')]/r0 * 1e6 * 1e-3
cRe0 = c_err[('0',0,'R')]/r0 * 1e6 * 1e-3
cL0 = c[('0',0,'L')]/r0 * 1e6 * 1e-3
cLe0 = c_err[('0',0,'L')]/r0 * 1e6 * 1e-3

eta_rw = 100 * (cL*r1) / (cL0*r0) # must multiply by ROI width again to get correct efficiency
eta_rw_err = eta_rw * np.sqrt((cLe/cL)**2 + (cLe0/cL0)**2)

N_in = (c[('0',0,'R')] + c[('0',0,'L')])*10
                          
print "Retrieval: %.2f +- %.2f kHz R %.2f +- %.2f kHz L\nInput: %.2f +- %.2f kHz R %.2f +- %.2f kHz L" % (cR, cRe, cL, cLe, cR0, cRe0, cL0, cLe0)

ext = (cR - dark_R)/(cL- dark_L)
ext_err = np.sqrt( (cRe/(cL-dark_L))**2 + (dark_R_err/(cL-dark_L))**2 + (ext*cLe/(cL-dark_L))**2 + (ext*dark_L_err/(cL - dark_L))**2 )
ext0 = (cR0 - dark_R)/(cL0- dark_L)
ext0_err = np.sqrt( (cRe0/(cL0-dark_L))**2 + (dark_R_err/(cL0-dark_L))**2 + (ext0*cLe0/(cL0-dark_L))**2 + (ext0*dark_L_err/(cL0 - dark_L))**2 )

ax = plt.subplot(223, axisbg=bgc)
ax.axvspan(*ROIs[0], color='white')
plt.errorbar(h_in[('L')].index.values, h_in[('L')], yerr=h_in_err[('L')].values, marker = 'o', color='#85c1fe', label="L")
plt.errorbar(h_in[('R')].index.values, h_in[('R')], yerr=h_in_err[('R')].values, marker = 'o', color='orange', label="R")
plt.ylabel(u"Zählrate (kHz)")
plt.xlabel(u"Zeit (\si{\micro\second})")
#plt.legend(title='Polarisation', fancybox=True)
#plt.axvline(x=ROIs[0][0], color='k', ls='solid')
#plt.axvline(x=ROIs[0][1], color='k', ls='solid')
ax.locator_params(tight=True, nbins=6)
plt.margins(0,.1)
plt.ylim(ymin=0)
#plt.text(10.6,175, r'$N_{\mathrm{in}} = \SI{%.1f}{}$' % N_in + '\n' + 
#                  r'$\epsilon_{L\mathrm{, in}} = \SI{%s}{\percent}$' % un2str(ext0*100, ext0_err*100), ha='center' )

h_L=pd.DataFrame(index=h_in[('L')].index.values-11.0,data=watt_factor*h_in[('L')].values)
h_L_err=pd.DataFrame(index=h_in[('L')].index.values-11.0,data=watt_factor*h_in_err[('L')].values)
h_R=pd.DataFrame(index=h_in[('R')].index.values-11.0,data=watt_factor*h_in[('R')].values)
h_R_err=pd.DataFrame(index=h_in[('R')].index.values-11.0,data=watt_factor*h_in_err[('R')].values)

plot_dict['223']={
    'A':{'type':'errorbar','y':h_L[0].to_json(),'yerr':h_L_err[0].to_json(),'xlabel':r'Zeit ($\mu s$)','ylabel':r'Leistung (fW)','num':'c','label':'L','xlim':(-0.6,0.2001),'ylim':(-watt_factor*30,800)} ,
    'B':{'type':'errorbar','y':h_R[0].to_json(),'yerr':h_R_err[0].to_json(),'label':'R'}
                                   
}


ax = plt.subplot(224, axisbg=bgc)
ax.axvspan(*ROIs[1], color='white')
plt.errorbar(h_ret[('L')].index.values, h_ret[('L')], yerr=h_ret_err[('L')].values, marker = 'o', color='#85c1fe', label='L')
plt.errorbar(h_ret[('R')].index.values, h_ret[('R')], yerr=h_ret_err[('R')].values, marker = 'o', color='orange', label='R')
#plt.ylabel(u"Zählrate (kHz)")
plt.xlabel(u"Zeit (\si{\micro\second})")
ax.locator_params(tight=True, nbins=6)
#plt.axvline(x=ROIs[1][0], color='k', ls='solid')
#plt.axvline(x=ROIs[1][1], color='k', ls='solid')
plt.margins(0,.1)
plt.ylim(ymin=0)
#plt.legend(title='Polarisation', fancybox=True)
#plt.text(15.65,10, r'$\epsilon_L = \SI{%s}{\percent}$' % (un2str(ext*100, ext_err*100)) + '\n' +
#                  r'$\eta_{\mathrm{rw}} = \SI{%s}{\percent}$' % un2str(eta_rw, eta_rw_err), ha='center' )

h_L=pd.DataFrame(index=h_ret[('L')].index.values-11.0,data=watt_factor*h_ret[('L')].values)
h_L_err=pd.DataFrame(index=h_ret[('L')].index.values-11.0,data=watt_factor*h_ret_err[('L')].values)
h_R=pd.DataFrame(index=h_ret[('R')].index.values-11.0,data=watt_factor*h_ret[('R')].values)
h_R_err=pd.DataFrame(index=h_ret[('R')].index.values-11.0,data=watt_factor*h_ret_err[('R')].values)

plot_dict['224']={
    'A':{'type':'errorbar','y':h_L[0].to_json(),'yerr':h_L_err[0].to_json(),'xlabel':r'Zeit ($\mu s$)','xlim':(4.0,5.5),'ylim':(-watt_factor*2,watt_factor*20),'num':'d','label':'L'} ,
    'B':{'type':'errorbar','y':h_R[0].to_json(),'yerr':h_R_err[0].to_json(),'label':'R'}
                                   
}

plt.tight_layout()
#savefig("extinction")
plt.show()

with io.open('memory_extinction.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))