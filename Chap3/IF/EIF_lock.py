# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 12:15:27 2016

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

def open_data(datei):
    print datei
    data=np.loadtxt(datei,skiprows=5, delimiter=',', converters = {0: lambda s: float(s.strip('"')), 1: lambda s: float(s.strip('"'))})
    return data.T

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
    
    
data_dir = 'InterferometerLockCharacterization'
lock = open_data(os.path.join(data_dir,'C1InLoop200ms00000.txt'))
no_lock = open_data(os.path.join(data_dir,'C1OutOfLoop00000.txt'))

f=plt.figure(0)
plt.subplot(121)
plt.xlabel('Zeit (s)')
plt.ylabel(r'cos($\Delta \phi$)')
plt.plot(no_lock[0], no_lock[1], 'k', rasterized=True)
plt.plot(lock[0], lock[1], 'g', rasterized=True)
plt.xlim(-1,1)
plt.ylim(-1.1, 1.1)

ax = plt.subplot(122)
ax.locator_params(tight=True, nbins=6)
plt.plot(times, results[4], 'r-')
plt.xlim(0,10)
plt.xlabel('Zeit (h)')
plt.ylabel('Azimut (rad)')

plt.tight_layout()
plt.show()

plot_dict={}

h2=pd.DataFrame(index=no_lock[0],data=no_lock[1])
h=pd.DataFrame(index=lock[0],data=lock[1])
plot_dict['121']={
    'A':{'type':'plot','y':h[0].to_json(),'xlabel':u'Zeit (s)','ylim':(-1.1,1.45),'xlim':(-1.,1.),'num':'a','ylabel':u'cos($\Delta \phi$)'},
    'B':{'type':'plot','y':h2[0].to_json()}
}


h=pd.DataFrame(index=times,data= results[4])
plot_dict['122']={
    'A':{'type':'plot','y':h[0].to_json(),'xlabel':u'Zeit (h)','ylim':(-0.47,-0.38),'xlim':(0,10.),'num':'b','ylabel':u'Azimut (rad)'}
}

with io.open('eif_lock.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))