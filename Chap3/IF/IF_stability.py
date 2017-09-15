# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 16:39:16 2016

@author: tstolz
"""

import numpy as np
import pylab as plt
import os
import pandas as pd
import json
import io

import quPy as qp
import Tomography as tom

data_dir = './stability_001/'

files = [data_dir + 'photons%d' % i for i in range(380,2200) ]

fields = [qp.MillisecondsSinceStart]
categorizers = [qp.POL]
from_header = ['TomographyBases', 'CycleID', 'LoadCycles']
profile = [fields, categorizers, {}, 0, 1000, from_header]
T = tom.Tomography(files[0:1], profile, [[]])

starttime = int(T.get_categories()['CycleID'].pop())
results = []
times = []
num_clicks = []
i=0
while i < len(files):
    print "file ", i
    T.read_files(files[i:i+3])
    i+=3
    while not T.get_categories()['TomographyBases'] == set(["22","11","00"]) and i < len(files):
        i+=1
        T.read_files(files[i:i+1])
    num_clicks.append( int(T.getCounts(level=[])[0]) )
    current_time = (np.mean(map(int, list(T.get_categories()['CycleID']))) - starttime)/3600.
    if current_time > 10:
        break
    results.append(T.getVectorCounts(level=['pol']))
    times.append(current_time)
    T.clear_data()

results=np.transpose(results)
        
print "Azimuth +- standarddeviation: %1.3f +- %1.3f" % ( np.mean(results[4]), np.std(results[4]) )
print "Polar angle +- stddev: %1.3f +- %1.3f" % ( np.mean(results[2]), np.std(results[2]) )
print "DOP +- stddev: %1.1f +- %1.1f %%" % ( 100*np.mean(results[0]), 100*np.std(results[0]) )

f=plt.figure(0)
ax = plt.subplot(311)
ax.locator_params(tight=True, nbins=6)
plt.plot(times, results[0], 'g-')
plt.xlim(0,10)
plt.ylabel('Polarisationsgrad')
ax.set_xticklabels([])

ax = plt.subplot(312)
ax.locator_params(tight=True, nbins=6)
plt.plot(times, results[2], 'b-')
plt.xlim(0,10)
plt.ylabel('Polarwinkel (rad)')
ax.set_xticklabels([])

ax = plt.subplot(313)
ax.locator_params(tight=True, nbins=6)
plt.plot(times, results[4], 'r-')
plt.xlim(0,10)
plt.xlabel('Zeit (h)')
plt.ylabel('Azimut (rad)')

plt.tight_layout()
plt.show()

