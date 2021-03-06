import locale
import numpy as np
import io
import json
import pandas as pd
import ast
import os
# Set to German locale to get comma decimal separater
locale.setlocale(locale.LC_NUMERIC, 'deu_deu')
#locale.setlocale(locale.LC_NUMERIC, 'de_DE.utf8')
import matplotlib as mpl
mpl.use('pgf')

preamble = [
    # use utf8 fonts becasue your computer can handle it :)
    r"\usepackage[utf8x]{inputenc}",
    # plots will be generated using this preamble
    r"\usepackage[T1]{fontenc}",
    r"\usepackage[ngerman]{babel}",
    r"\usepackage{siunitx}",
    r"\usepackage{lmodern}",
    r"\usepackage{amsmath}",
    r"\usepackage{amsfonts}",
    r"\sisetup{detect-all}",
    r"\sisetup{locale = DE}"
]


pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "text.latex.unicode": True,
    "font.family": "sans-serif",
    # blank entries should cause plots to inherit fonts from the document
    "font.serif": [],
    "font.sans-serif": ['Helvetica'],
    "font.monospace": [],
    "axes.labelsize": 11,               # LaTeX default is 10pt font.
    "font.size": 11,
    "legend.fontsize": 10,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.figsize": [0.9*5.67, 1.76],   # default fig size of 0.9 textwidth
    "errorbar.capsize": 0,             # set standard
    "markers.fillstyle": 'none',
    "lines.markersize": 1,
    "lines.linewidth": 1.5,
    "legend.fancybox": True,
    "mathtext.fontset": "cm",
    "text.latex.preamble": preamble,
    #    "pgf.debug" : True,
    #"legend.loc": 1,
    "pgf.preamble": preamble,
    "legend.numpoints": 1,
    "legend.scatterpoints": 1,
    "axes.formatter.use_locale": True,
    "figure.subplot.bottom" : 0.19
}

mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt


class ThesisPlot(object):
    
    colors=['b','r','g','k','y']
    linestyles=['-','-','-','-','-','-']
    linewidths=[1.5,1.5,1.5,1.5,1.5,1.5]
    markers=['o','o','o','o','o']
    elinewidth = 0.8

    def __init__(self):
        self.dicts=dict()

    def generatePlots(self):
        for d in self.dicts:
            self.f=plt.figure()
            
#            self.dicts[d]['json']=list(self.dicts[d]['json'])[::-1]
#            print type(self.dicts[d]['json'])
            
            for subplot_i, sp_ax in zip(range(len(self.dicts[d]['json'])-1,-1,-1), self.dicts[d]['json']):
                ax=self.f.add_subplot(sp_ax)
                
                ylabeloffset=self.dicts[d]['ylabeloff'][subplot_i]
                print "ylavel off: "
                print ylabeloffset
                
                if ylabeloffset:
                    ax.yaxis.set_label_coords(ylabeloffset,0.5)
                
                num_curves = len(self.dicts[d]['json'][sp_ax])
                
                cs=self.dicts[d]['color']
                if cs==None:
                    cs=self.colors
                if len(np.shape(cs))==2:
                    cs=cs[subplot_i]
                    
                
                ls=self.dicts[d]['linestyle']
                if ls==None:
                    ls=self.linestyles
                if len(np.shape(ls))==2:
                    ls=ls[subplot_i]
                    
                lws=self.dicts[d]['linewidths']
                if lws==None:
                    lws=self.linewidths
                if len(np.shape(lws))==2:
                    lws=lws[subplot_i]
                    
                tl=self.dicts[d]['tight']
                if tl:
                    self.f.tight_layout(w_pad=self.dicts[d]['wpad'],h_pad=self.dicts[d]['hpad'])
                    
                xticks = self.dicts[d]['xticks']
                print "shape"
                print np.shape(xticks)
                if np.shape(xticks)[0]==1:
                    ax.xaxis.set_ticks(xticks[0])
                elif len(np.shape(xticks))==1 and len(xticks) != 0:
                    if xticks[subplot_i] is not None:
                        ax.xaxis.set_ticks(xticks[subplot_i])
                
                    
                yticks = self.dicts[d]['yticks']
                if len(np.shape(yticks))==0:
                    ax.yaxis.set_ticks(yticks)
                elif len(np.shape(yticks))==1 and  len(yticks) != 0:
                    if yticks[subplot_i] is not None:
                        ax.yaxis.set_ticks(yticks[subplot_i])
                
                    
                ms=self.dicts[d]['markers']
                if ms==None:
                    ms=self.markers
                if len(np.shape(ms))==2:
                    ms=ms[subplot_i]
                
                
                
                for (c,l,lw,sp,m) in zip(cs[:num_curves],ls[:num_curves],lws[:num_curves],self.dicts[d]['json'][sp_ax],ms[:num_curves]):
                    
                    if self.dicts[d]['json'][sp_ax][sp]['type']=='errorbar':
                        df = pd.DataFrame.from_dict(json.loads(self.dicts[d]['json'][sp_ax][sp]['y']),orient='index')
                        df.set_index(np.array(df.index.values,dtype=np.float32),inplace=True)
                        df.sort(inplace=True)
                        x=np.array(df.index.values,dtype=np.float32)
                        y=np.array(df.values,dtype=np.float32)
                        
                        dfErr = pd.DataFrame.from_dict(json.loads(self.dicts[d]['json'][sp_ax][sp]['yerr']),orient='index')
                        dfErr.set_index(np.array(dfErr.index.values,dtype=np.float32),inplace=True)
                        dfErr.sort(inplace=True)
                        yerr=np.array(dfErr.values,dtype=np.float32)
                        
                        try:
                            label=self.dicts[d]['json'][sp_ax][sp]['label']
                        except:
                            label=None
                            
                        try:
                            ax.set_xlabel(self.dicts[d]['json'][sp_ax][sp]['xlabel'])
                        except:
                            print "No xlabel in %s ax %s"%(d, sp_ax)
                            
                        try:
                            ax.set_ylabel(self.dicts[d]['json'][sp_ax][sp]['ylabel'])
                        except:
                            print "No xlabel in %s ax %s"%(d, sp_ax)
                            
                        try:
                            xl=self.dicts[d]['json'][sp_ax][sp]['xlim']
                            ax.set_xlim(*xl)
                        except:
                            print "No x-limit found"
                            
                        try:
                            yl=self.dicts[d]['json'][sp_ax][sp]['ylim']
                            ax.set_ylim(*yl)
                        except:
                            print "No y-limit found"
                            
                        
                        ax.errorbar(x,y,yerr=yerr,label=label,color=c,ls=l,lw=lw,marker=m,markersize='5', elinewidth=self.elinewidth)
                    elif self.dicts[d]['json'][sp_ax][sp]['type']=='plot':
                        df = pd.DataFrame.from_dict(json.loads(self.dicts[d]['json'][sp_ax][sp]['y']),orient='index')
                        df.set_index(np.array(df.index.values,dtype=np.float32),inplace=True)
                        df.sort(inplace=True)
                        x=np.array(df.index.values,dtype=np.float32)
                        y=np.array(df.values,dtype=np.float32)
                        
                        try:
                            m=self.dicts[d]['json'][sp_ax][sp]['margin']
                            ax.margins(*m)
                            print "found margin"
                        except:
                            print "No margin"
                        
                        try:
                            label=self.dicts[d]['json'][sp_ax][sp]['label']
                        except:
                            label=None
                            
                        try:
                            ax.set_xlabel(self.dicts[d]['json'][sp_ax][sp]['xlabel'])
                        except:
                            print "No xlabel in %s ax %s"%(d, sp_ax)
                            
                        try:
                            ax.set_ylabel(self.dicts[d]['json'][sp_ax][sp]['ylabel'])
                        except:
                            print "No xlabel in %s ax %s"%(d, sp_ax)
                         
                        try:
                            xl=self.dicts[d]['json'][sp_ax][sp]['xlim']
                            ax.set_xlim(*xl)
                        except:
                            print "No x-limit found"
                        
                        try:
                            yl=self.dicts[d]['json'][sp_ax][sp]['ylim']
                            ax.set_ylim(*yl)
                        except:
                            print "No y-limit found"
                        
                        ax.plot(x,y,label=label,color=c,ls=l,lw=lw)
                    elif self.dicts[d]['json'][sp_ax][sp]['type']=='scatter':
                        df = pd.DataFrame.from_dict(json.loads(self.dicts[d]['json'][sp_ax][sp]['y']),orient='index')
                        df.set_index(np.array(df.index.values,dtype=np.float32),inplace=True)
                        df.sort(inplace=True)
                        x=np.array(df.index.values,dtype=np.float32)
                        y=np.array(df.values,dtype=np.float32)
                        
                        try:
                            m=self.dicts[d]['json'][sp_ax][sp]['margin']
                            ax.margins(*m)
                            print "found margin"
                        except:
                            print "No margin"
                        
                        try:
                            label=self.dicts[d]['json'][sp_ax][sp]['label']
                        except:
                            label=None
                            
                        try:
                            ax.set_xlabel(self.dicts[d]['json'][sp_ax][sp]['xlabel'])
                        except:
                            print "No xlabel in %s ax %s"%(d, sp_ax)
                            
                        try:
                            ax.set_ylabel(self.dicts[d]['json'][sp_ax][sp]['ylabel'])
                        except:
                            print "No xlabel in %s ax %s"%(d, sp_ax)
                         
                        try:
                            xl=self.dicts[d]['json'][sp_ax][sp]['xlim']
                            ax.set_xlim(*xl)
                        except:
                            print "No x-limit found"
                        
                        try:
                            yl=self.dicts[d]['json'][sp_ax][sp]['ylim']
                            ax.set_ylim(*yl)
                        except:
                            print "No y-limit found"
                        
                        ax.scatter(x,y,label=label,color=c,marker='o',s=16)
                    elif self.dicts[d]['json'][sp_ax][sp]['type']=='axh':
                        try:
                            label=self.dicts[d]['json'][sp_ax][sp]['label']
                        except:
                            label=None
                            
                        ax.axhline(np.float(self.dicts[d]['json'][sp_ax][sp]['y']),label=label,color=c,ls=l,lw=lw)
                    elif self.dicts[d]['json'][sp_ax][sp]['type']=='axv':
                        try:
                            label=self.dicts[d]['json'][sp_ax][sp]['label']
                        except:
                            label=None
                        ax.axvline(np.float(self.dicts[d]['json'][sp_ax][sp]['y']),label=label,color=c,ls=l,lw=lw)
                        
#                    
                    if self.dicts[d]['legend']:
                        ax.legend(loc=self.dicts[d]['loc'])
                    try:
                        num=self.dicts[d]['json'][sp_ax][sp]['num']
                    except:
                        num=None
                    if num is not None:
                        ax.text(0.1, 0.9, r'\textbf{(' + num + ')}', transform=ax.transAxes,
                                weight='bold', ha='center', va='center')
                        xwin, ywin = ax.transAxes.transform((0.1, 0.9))
                        for l in ax.yaxis.get_major_ticks():
                            # check if a label overlaps with enumeration
                            bbox = l.label1.get_window_extent()
                            print bbox, xwin, ywin
                            if self._overlaps(np.array(bbox), xwin, ywin):
                                l.label1.set_visible(False)
                    
#            if len(self.f.axes) > 1:
#                for n, ax in enumerate(self.f.axes):
#                    ax.text(0.1, 0.9, r'\textbf{(' + chr(len(self.f.axes)-1-n + 97) + ')}', transform=ax.transAxes,
#                            weight='bold', ha='center', va='center')
#                    # label position in window coordinates
#                    xwin, ywin = ax.transAxes.transform((0.1, 0.9))
#                    for l in ax.yaxis.get_major_ticks():
#                        # check if a label overlaps with enumeration
#                        bbox = l.label1.get_window_extent()
#                        print bbox, xwin, ywin
#                        if self._overlaps(np.array(bbox), xwin, ywin):
#                            l.label1.set_visible(False)

            
            s=self.figsize(self.dicts[d]['size'],1.0)
            self.f.subplots_adjust(bottom=self.dicts[d]['bottom']) 
            self.f.set_size_inches(*s)
            print self.dicts[d]['outfile']
            self.f.savefig(self.dicts[d]['outfile'])
            self.f.savefig(self.dicts[d]['outfile']+".pdf")
#            self.f.clear()
    
    def _overlaps(self, bbox, x, y, dist=10):
        xs, ys = bbox.T
        if (np.min(np.abs(xs - x)) > dist and np.prod(xs - x) > 0) or \
                (np.min(np.abs(ys - y)) > dist and np.prod(ys - y) > 0):
            return False
        else:
            # print np.min(np.abs(xs-x)), np.prod(xs-x), np.min(np.abs(ys-y)),
            # np.prod(ys-y)
            return True
        
    def parsePlotDict(self,filename):
        with io.open(filename, 'r', encoding='utf-8') as f:
            plotDict=json.load(f)
            
        return plotDict

    def addPlot(self,name,outname,figid,size=2,ls=None,cs=None,lw=None,tl=False,w_pad=2.,h_pad=2.,legend=False,lloc=1,m=None, xticks=[], yticks=[], bottom=0.2,yoffset=None):
        nplots=len(self.parsePlotDict(name))
        if yoffset is None:
            yoffset=[None for _ in range(nplots)]
                    
        self.dicts.update({figid:{'infile':name,
                                  'outfile':outname,
                                  'size':size,
                                  'json':self.parsePlotDict(name),
                                  'linestyle':ls,
                                  'color':cs,
                                  'linewidths':lw,
                                  'markers':m,
                                  'tight':tl,
                                  'wpad':w_pad,
                                  'hpad':h_pad,
                                  'loc':lloc,
                                  'legend':legend,
                                  'xticks':xticks,
                                  'yticks':yticks,
                                  'bottom':bottom,
                                  'ylabeloff':yoffset}})

    def figsize(self, rows, scale):
        # Get this from LaTeX using \the\textwidth
        fig_width_pt = 405.45183
        inches_per_pt = 1.0 / 72.27                       # Convert pt to inch
        # Aesthetic ratio (you could change this) * 0.5
        golden_mean = rows * 0.51 * (np.sqrt(5.0) - 1.0) / 2.0
        fig_width = fig_width_pt * inches_per_pt * scale    # width in inches
        fig_height = fig_width * golden_mean              # height in inches
        fig_size = [fig_width, fig_height]
        return fig_size

if __name__=='__main__':
    TP=ThesisPlot()
#    TP.addPlot("Chap2\Groupdelay\groupdelay.json","2_2_groupdelay.pgf","Chap2_Fig2.2")
#    TP.addPlot(os.path.join("Chap2","Suszept","suszept.json"),"2_1_suszept.pgf","Chap2_Fig2.1",size=2,cs=['b','k','r'],ls=['-','--','-'],lw=[1.5,1,1.5], bottom=0.11)
#    TP.addPlot(r"blank1.json","2_3_blank.pgf","Chap2_Fig2.3",size=1)
    TP.addPlot(os.path.join("Chap2","Transient","eit_propagation.json"),"2_3_eit_propagation.pgf","Chap2_Fig2.3",size=1,tl=True,yoffset=[-0.13,None],w_pad=1.6, bottom=0.22)
#    TP.addPlot("Chap2\Foerster\defect.json","2_4_foerster_defect.pgf","Chap2_Fig2.4",size=1.0,legend=True,cs=['b','r','k'])
#    TP.addPlot(os.path.join("Chap2","BlockadeSuszept","blockade_suszept.json"),"2_7_blockade_suszept.pgf","Chap2_Fig2.7",size=1.0,legend=True,cs=['b','r'], bottom=0.22,h_pad=0.0,w_pad=1.2,tl=True)
#    TP.addPlot(os.path.join("Chap2","KondPhase","cond_phase.json"),"2_8_cond_phase.pgf","Chap2_Fig2.8",size=1.0,cs=['b','k','r'],yoffset=[-0.145,None],legend=True,h_pad=0.0,w_pad=1.4,tl=True,lloc=9, bottom=0.22)
#    TP.addPlot(os.path.join("Chap2","MoleculeMemory","memory.json"),"2_10_memory.pgf","Chap2_Fig2.10",size=1.0,cs=['b','r'],yoffset=[-0.19,None],legend=True,h_pad=0.0,w_pad=1.9,tl=True,lloc=4, bottom=0.22)
#    TP.addPlot(os.path.join("Chap3","Plugs","density_profile.json"),"3_2_density_profile.pgf","Chap3_Fig3.2",size=1.0,cs=['b','r'],legend=False,h_pad=0.0,w_pad=1.5,tl=False,lloc=4, bottom=0.22)
#    TP.addPlot(os.path.join("Chap3","Plugs","pluglength.json"),"3_1_pluglength.pgf","Chap3_Fig1.1",size=1.0,cs=['b','r'],legend=False,h_pad=0.0,w_pad=1.5,tl=True,lloc=4, bottom=0.22)
#    TP.addPlot(os.path.join("Chap3","IF","eif_lock.json"),"3_10_eif_lock.pgf","Chap3_Fig3.10",size=1.0,cs=['b','r'],legend=True,h_pad=0.0,w_pad=2.3,tl=True,lloc=4, bottom=0.22)
#    TP.addPlot(os.path.join("Chap3","Laser","cavity_characterization.json"),"3_5_cavity.pgf","Chap3_Fig3.5",size=1.0,cs=['b','r'],legend=True,h_pad=0.0,w_pad=1.,tl=True,lloc=1, bottom=0.22)
#    TP.addPlot(r"Chap5\sideband_postselected_phaseshift.json","5_2_phaseshift.pgf","Chap5_Fig5.2",size=1.0,cs=['b','r'],legend=True,h_pad=0.0,w_pad=1.,tl=True,lloc=1)
#    TP.addPlot(os.path.join("Chap5","sideband_postselected_phaseshift.json"),"5_2_phaseshift.pgf","Chap5_Fig5.2",size=1.0,cs=['b','k','r'],legend=True,h_pad=0,w_pad=0.,lloc=1, bottom=0.22)
#    TP.addPlot(os.path.join("Chap5","spectrum.json"),"5_1_spectrum.pgf","Chap5_Fig5.2",size=1.0,cs=['b','r','r','b'], yoffset=[-0.18,None],legend=True,h_pad=0,w_pad=2.,lloc=1,tl=True, bottom=0.22)
#    TP.addPlot(os.path.join("Chap5","pol_spectra.json"),"5_4_spectra.pgf","Chap5_Fig5.4",size=1.0,cs=['b','b','r','k','r'],yoffset=[-0.14,None],legend=False,tl=True,h_pad=0,w_pad=1.8,lloc=1,ls=['-','','-','-',''],m=['','o','','','o'], bottom=0.22)
#    TP.addPlot(os.path.join("Chap5","cond_phase_vs_density.json"),"5_7_phase.pgf","Chap5_Fig5.7",size=1.0,cs=['b','b','r','r'],legend=False,tl=False,h_pad=0,w_pad=1.8,lloc=1,ls=['-','','-',''],m=['','o','','o',''], bottom=0.23)
#    TP.addPlot(os.path.join("Chap5","propagation.json"),"5_5_propagation.pgf","Chap5_Fig5.5",size=1.0,cs=['b','b','r','r'],xticks=[np.arange(11,14.5,1)],yoffset=[-0.18,None],legend=False,tl=True,h_pad=0,w_pad=2.2,lloc=1,ls=['-','','-',''],m=['','o','','o',''], bottom=0.22)
#    TP.addPlot(os.path.join("Chap5","storage_retrieval.json"),"5_6_storage.pgf","Chap5_Fig5.6",size=1.0,cs=['b','r'],legend=False,tl=False,h_pad=0,w_pad=0,lloc=1,ls=['','-'], bottom=0.22)
#    TP.addPlot(os.path.join("Chap2","Molecules","avg_number.json"),"2_2_moleculetest.pgf","Chap5_Fig2.2",size=1.0)
#    TP.addPlot(os.path.join("Chap6","memory_spectra.json"),"6_1_spectra.pgf","Chap6_Fig6.1",size=2,cs=[['b','b'],['r','r'],['b','b'],['r','r']],ls=['-',''], tl=True,h_pad=-1.,w_pad=0,yticks=[None,None,[-9,-6,-3,0,3,6,9],None], bottom=0.11)
#    TP.addPlot(os.path.join("Chap6","memory_extinction.json"),"6_2_extinction.pgf","Chap6_Fig6.2",size=2,cs=['b','r'],xticks=[np.arange(-0.6,0.4,0.2),None,np.arange(-0.6,0.4,0.2),None],yticks=[None,None,np.arange(0,900,200),None], bottom=0.11)
#    TP.addPlot(os.path.join("Chap6","memory_coherence.json"),"6_3_coherence.pgf","Chap6_Fig6.3",size=3,cs=['b','b','r','r'], ls=['','-','','-'], bottom=0.075)
#    TP.addPlot(os.path.join("Chap6","darktime.json"),"6_4_darktime.pgf","Chap6_Fig6.4",size=1,yoffset=[-0.14,None],cs=(('b','b','b'),('r','b','b')), ls=[['','-',''],['','-','']], tl=True,h_pad=0,w_pad=1.5, bottom=0.22)
#    TP.addPlot(os.path.join("Chap2","Molecules","avg_number.json"),"2_avgnumber.pgf","Chap2_Fig2.2")    
    TP.generatePlots()