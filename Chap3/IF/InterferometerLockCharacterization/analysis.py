# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 12:15:27 2016

@author: tstolz
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import os
mpl.rcParams['axes.color_cycle']=['k', 'b', 'r', 'g', 'c', 'm', 'y']
mpl.rcParams["savefig.directory"] = os.getcwd()
import numpy as np
from scipy.optimize import curve_fit
import Tkinter as tk
import tkFileDialog as tkf

def sine(t,o,a,f,p):
    return o+a*np.sin(2*np.pi*f*t+p)

def fourier(data):
    time, volt = data    
    N=len(time)
    dt=(time[-1]-time[0])/(float(N)-1)
    fft = np.abs(np.fft.rfft(volt))
    freqaxis= np.fft.rfftfreq(N,dt)
#    plt.figure()
    plt.plot(freqaxis, fft)
    plt.xlabel("frequency (Hz)")
    plt.ylabel("amplitude (a.u.)")
    plt.show()
    return freqaxis, fft

def fit_sine(data, freq=None):
    ampl = (data[1].max()-data[1].min())/2
    offs = data[1].min() + ampl
    time, volt = data    
    if not freq:
        N=len(time)
        dt=(time[-1]-time[0])/(float(N)-1)
        fft = np.abs(np.fft.rfft(volt))
        freqaxis= np.fft.rfftfreq(N,dt)
        freq=freqaxis[5 + fft[5:].argmax()]
        plt.figure()
        plt.plot(freqaxis, fft)
        plt.xlabel("frequency (Hz)")
        plt.ylabel("amplitude (a.u.)")
        plt.show()
    phase=np.arcsin((data[1][0]-offs)/ampl)+2*np.pi*freq*data[0][0]
    if np.max(data[1]-sine(data[0], offs, ampl, freq, phase))>(0.9*ampl):
        phase=np.pi-phase
    plt.figure()
    plt.plot(data[0], data[1], data[0], sine(data[0], offs, ampl, freq, phase))
    plt.title("guessed parameters")
    print "guessed parameters: offs=%f, ampl=%f, freq=%f, phase=%f" % (offs, ampl, freq, phase)
    popt, pcov = curve_fit( sine , time, volt, [offs, ampl, freq, phase])
    print "fit result: offs=%f, ampl=%f, freq=%f, phase=%f" % tuple(popt)
    #popt = [offs,ampl,freq,phase]
    plt.figure()
    plt.xlabel('time (s)')
    plt.ylabel('amplitude (V)')
    plt.grid()
    plt.plot(data[0], data[1], label='data')
    plt.plot(data[0], sine(data[0], *popt), label='sine fit')
    plt.legend()
    return popt, np.sqrt(np.diag(pcov))
    
def open_data(datei):
    print datei
    data=np.loadtxt(datei,skiprows=5, delimiter=',', converters = {0: lambda s: float(s.strip('"')), 1: lambda s: float(s.strip('"'))})
    return data.T
    
def get_amplitude(data,win_size=100):
    '''analyzes the peak to peak range of the timetrace "data" within windows 
    of size "win_size". timestamps are calculated as mean values of the window,
    assuming evenly spaced sampling points.'''
    time=[]
    ampl=[]
    for i in xrange(0,len(data[0]),win_size):
        time.append(data[0,i:i+win_size].mean())
        ampl.append(data[1,i:i+win_size].max()-data[1,i:i+win_size].min())
    return np.array([time,ampl])

def chop(data, offs=0, interval=12500, clear=200):
    '''chops the dataset into two sets containing blocks of length "interval",
    skipping "clear" datapoints around the cuts and placing the first cut at "offs"'''
    time1, time2, ampl1, ampl2 = ([],[],[],[])
    for i in xrange(len(data[0])):
        if (i-offs)%interval > (interval-clear/2) or (i-offs)%interval < clear/2:
            continue
        elif ((i-offs)/interval) % 2 == 0:
            time1.append(data[0][i])
            ampl1.append(data[1][i])
        elif ((i-offs)/interval) % 2 == 1:
            time2.append(data[0][i])
            ampl2.append(data[1][i])
    return time1, time2, ampl1, ampl2

def phase_plot(ref, errsig):
    '''takes a reference signal (sinusoidal, providing phase information), 
    normalizes time and amplitude axes and performs an arcsin transformation 
    to obtain a plot of the phase evolution over time'''
    ref=np.asarray(ref)
    errsig=np.asarray(errsig)
    #normalize the time axes
    tr=(ref[0]-ref[0][0])/(ref[0].max()-ref[0].min())
    t=(errsig[0]-errsig[0][0])/(errsig[0].max()-errsig[0].min())
    #normalize and center the reference signal
    ar=2*(ref[1]-(ref[1].max()-ref[1].min())/2-ref[1].min())/(ref[1].max()-ref[1].min())
    #do the same with the error signal
    a=2*(errsig[1]-(ref[1].max()-ref[1].min())/2-ref[1].min())/(ref[1].max()-ref[1].min())
    #calculate phases    
    pr=np.arcsin(ar)
    p=np.arcsin(a)
    plt.figure()
    plt.xlabel('normalized time')
    plt.ylabel('phase (rad)')
    plt.plot(tr, pr, 'g.', t, p, 'r.')
    plt.show()
    return np.array([tr, ar, pr, t, a, p])

    
if __name__ == "__main__":
    
    plt.ion()
    lock = open_data('C1InLoop200ms00000.txt')
    no_lock = open_data('C1OutOfLoop00000.txt')

    plt.figure(figsize=(7,4))
    plt.xlabel('time (s)')
    plt.ylabel('cos($\phi$)')
    plt.grid()
    plt.plot(no_lock[0], no_lock[1])
    plt.plot(lock[0], lock[1], 'b')
    plt.xlim(-1,1)
    plt.ylim(-1.1, 1.1)
    plt.tight_layout()
    
    plt.figure(figsize=(7,4))
    plt.xlabel('time (s)')
    plt.ylabel('cos($\phi$)')
    plt.grid()
    plt.plot(no_lock[0], no_lock[1])
    plt.xlim(-1,1)
    plt.ylim(-1.1, 1.1)
    plt.tight_layout()
    