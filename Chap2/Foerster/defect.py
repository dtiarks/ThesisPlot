# -*- coding: utf-8 -*-
"""
Created on Fri Jun 02 16:27:16 2017

@author: daniel
"""

import Tomography as tom
import quPy as qp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import json
import io


dataNN=np.loadtxt("foersterdefect_n_n.tsv")
dataNN_2=np.loadtxt("foersterdefect_n_n-2.tsv")

plot_dict={}

f=plt.figure()
#plt.subplot(121)
#plt.plot(dataNN[:,0], dataNN[:,1]/1e9, marker='o',ls='',color='b')
#plt.plot(dataNN[:,0], dataNN[:,2]/1e9, marker='o',ls='',color='r')
#plt.plot(dataNN[:,0], dataNN[:,3]/1e9, marker='o',ls='',color='g')
#plt.plot(dataNN[:,0], dataNN[:,4]/1e9, marker='o',ls='',color='orange')
#plt.xlim((45,120))
#plt.ylim((-4,0))


h=pd.DataFrame(index=dataNN[:,0],data=dataNN[:,1]/1e9)
h2=pd.DataFrame(index=dataNN[:,0],data=dataNN[:,2]/1e9)
h3=pd.DataFrame(index=dataNN[:,0],data=dataNN[:,3]/1e9)
h4=pd.DataFrame(index=dataNN[:,0],data=dataNN[:,4]/1e9)
#plot_dict['121']={
#    'A':{'type':'scatter','y':h[0].to_json(),'label':u'$nP_{3/2}(n-1)P_{3/2}$','ylabel':u'$\Delta E_p/h$ (GHz)','xlabel':u'Hauptquantenzahl $n$','xlim':(45,120),'ylim':(-4,0)},                
#    'B':{'type':'scatter','y':h2[0].to_json(),'label':u'$nP_{1/2}(n-1)P_{3/2}$'},
#    'C':{'type':'scatter','y':h3[0].to_json(),'label':u'$nP_{3/2}(n-1)P_{1/2}$'},
#    'D':{'type':'scatter','y':h4[0].to_json(),'label':u'$nP_{1/2}(n-1)P_{1/2}$'}
#}

plt.subplot(111)
plt.plot(dataNN_2[:,0], dataNN_2[:,1]/1e6, marker='o',ls='',color='b')
plt.plot(dataNN_2[:,0], dataNN_2[:,2]/1e6, marker='o',ls='',color='r')
plt.plot(dataNN_2[:,0], dataNN_2[:,3]/1e6, marker='o',ls='',color='g')
plt.plot(dataNN_2[:,0], dataNN_2[:,4]/1e6, marker='o',ls='',color='orange')
plt.xlim((64.5,76.5))
plt.axhline(0,c='k',lw=1.5)
plt.ylim((-100,100))

h=pd.DataFrame(index=dataNN_2[:,0],data=dataNN_2[:,1]/1e9)
h2=pd.DataFrame(index=dataNN_2[:,0],data=dataNN_2[:,2]/1e9)
h3=pd.DataFrame(index=dataNN_2[:,0],data=dataNN_2[:,3]/1e9)
h4=pd.DataFrame(index=dataNN_2[:,0],data=dataNN_2[:,4]/1e9)
plot_dict['111']={
    'A':{'type':'scatter','y':h2[0].to_json(),'label':u'$(n-1)P_{1/2}(n-2)P_{3/2}$','xlabel':u'Hauptquantenzahl $n$','ylabel':u'$\Delta E_p/h$ (GHz)','xlim':(64.5,76.5),'ylim':(-0.1,0.1)},                
    'B':{'type':'scatter','y':h3[0].to_json(),'label':u'$(n-1)P_{3/2}(n-2)P_{1/2}$'},
#    'C':{'type':'scatter','y':h3[0].to_json(),'label':u'$(n-1)P_{3/2}(n-2)P_{1/2}$'},
#    'D':{'type':'scatter','y':h4[0].to_json(),'label':u'$(n-1)P_{1/2}(n-2)P_{1/2}$'},
    'E':{'type':'axh','y':0}
}


with io.open('defect.json', 'w+') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False,indent=4)))


plt.show()