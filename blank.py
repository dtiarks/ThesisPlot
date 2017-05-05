# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:13:15 2017

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

h=pd.DataFrame(index=[0],data=[0])

#f=plt.figure()
#
#ax1=f.add_subplot(221)
#ax1.plot([0],[0])
#
#ax1=f.add_subplot(222)
#ax1.plot([0],[0])
#
#ax1=f.add_subplot(223)
#ax1.plot([0],[0])
#
#ax1=f.add_subplot(224)
#ax1.plot([0],[0])

plot_dict={}
plot_dict['121']={
                    'A':{'type':'plot','label':'test','y':h.to_json()}
                 }
                
plot_dict['122']={
                    'A':{'type':'plot','label':'test','y':h.to_json()}
                 }
                 
#plot_dict['222']={
#                    'A':{'type':'plot','label':'test','y':h.to_json()}
#                 }
#                 
#plot_dict['223']={
#                    'A':{'type':'plot','label':'test','y':h.to_json()}
#                 }
#                 
#plot_dict['224']={
#                    'A':{'type':'plot','label':'test','y':h.to_json()}
#                 }
                 
with io.open('blank1.json', 'w') as f:
  f.write(unicode(json.dumps(plot_dict, ensure_ascii=False)))