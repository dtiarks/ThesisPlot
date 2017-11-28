import locale
import numpy as np
import io
import json
import pandas as pd
import ast
import os
# Set to German locale to get comma decimal separater
#locale.setlocale(locale.LC_NUMERIC, 'deu_deu')
locale.setlocale(locale.LC_NUMERIC, 'de_DE.utf8')
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
    "figure.figsize": [1*5.67, 1.76],   # default fig size of 0.9 textwidth
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

    def __init__(self):
        self.dicts=dict()

    def generatePlots(self):
        for d in self.dicts:
            self.f=plt.figure()
            
            for sp_ax in self.dicts[d]['json']:
                ax=self.f.add_subplot(sp_ax)
                
                cs=self.dicts[d]['color']
                if cs==None:
                    cs=self.colors
                
                ls=self.dicts[d]['linestyle']
                if ls==None:
                    ls=self.linestyles
                    
                lws=self.dicts[d]['linewidths']
                if lws==None:
                    lws=self.linewidths
                    
                tl=self.dicts[d]['tight']
                if tl:
                    self.f.tight_layout(w_pad=self.dicts[d]['wpad'],h_pad=self.dicts[d]['hpad'])
                    
                ms=self.dicts[d]['markers']
                print "markers: %s"%ms
                if ms==None:
                    ms=self.markers
                
                
                for (c,l,lw,sp,m) in zip(cs[:len(self.dicts[d]['json'][sp_ax])],ls[:len(self.dicts[d]['json'][sp_ax])],lws[:len(self.dicts[d]['json'][sp_ax])],self.dicts[d]['json'][sp_ax],ms[:len(self.dicts[d]['json'][sp_ax])]):
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
                            
                        print "marker: %s"%m
                        m='o'
                        
                        ax.errorbar(x,y,yerr=yerr,label=label,color=c,ls=l,lw=lw,marker=m,markersize='5')
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
            self.f.subplots_adjust(bottom=0.2) 
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

    def addPlot(self,name,outname,figid,size=2,ls=None,cs=None,lw=None,tl=False,w_pad=2.,h_pad=2.,legend=False,lloc=1,m=None):
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
                                  'legend':legend}})

    def figsize(self, rows, scale):
        # Get this from LaTeX using \the\textwidth
        fig_width_pt = 455.24416
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
#    TP.addPlot("Chap2\Suszept\suszept.json","2_1_suszept.pgf","Chap2_Fig2.1",size=2,cs=['b','k','r'],ls=['-','--','-'],lw=[1.5,1,1.5])
#    TP.addPlot(r"blank1.json","2_3_blank.pgf","Chap2_Fig2.3",size=1)
#    TP.addPlot("Chap2\Transient\eit_propagation.json","2_3_eit_propagation.pgf","Chap2_Fig2.3",size=1,tl=True,w_pad=1.4)
#    TP.addPlot("Chap2\Foerster\defect.json","2_4_foerster_defect.pgf","Chap2_Fig2.4",size=1.0,legend=True,cs=['b','r','k'])
#    TP.addPlot(r"Chap2\BlockadeSuszept\blockade_suszept.json","2_7_blockade_suszept.pgf","Chap2_Fig2.7",size=1.0,legend=True,cs=['b','r'])
#    TP.addPlot(r"Chap2\KondPhase\cond_phase.json","2_8_cond_phase.pgf","Chap2_Fig2.8",size=1.0,cs=['b','k','r'],legend=True,h_pad=0.0,w_pad=1.0,tl=True,lloc=9)
#    TP.addPlot(r"Chap2\Molecules\avg_number.json","2_9_avg_number.pgf","Chap2_Fig2.9",size=1.0,cs=['b','k','r'],legend=True,h_pad=0.0,w_pad=1.0,tl=True,lloc=2)
#    TP.addPlot(r"Chap2\MoleculeMemory\memory.json","2_10_memory.pgf","Chap2_Fig2.10",size=1.0,cs=['b','r'],legend=True,h_pad=0.0,w_pad=1.5,tl=True,lloc=4)
#    TP.addPlot(r"Chap3\Plugs\pluglength.json","3_1_pluglength.pgf","Chap3_Fig1.1",size=1.0,cs=['b','r'],legend=False,h_pad=0.0,w_pad=1.5,tl=True,lloc=4)
#    TP.addPlot(r"Chap3\Plugs\density_profile.json","3_2_density_profile.pgf","Chap3_Fig3.2",size=1.0,cs=['b','r'],legend=False,h_pad=0.0,w_pad=1.5,tl=False,lloc=4)
#    TP.addPlot(r"Chap3\IF\eif_lock.json","3_10_eif_lock.pgf","Chap3_Fig3.10",size=1.0,cs=['b','r'],legend=True,h_pad=0.0,w_pad=2.3,tl=True,lloc=4)
#    TP.addPlot(r"Chap3\Laser\cavity_characterization.json","3_5_cavity.pgf","Chap3_Fig3.5",size=1.0,cs=['b','r'],legend=True,h_pad=0.0,w_pad=1.,tl=True,lloc=1)
#    TP.addPlot(r"Chap5\sideband_postselected_phaseshift.json","5_2_phaseshift.pgf","Chap5_Fig5.2",size=1.0,cs=['b','r'],legend=True,h_pad=0.0,w_pad=1.,tl=True,lloc=1)
#    TP.addPlot(os.path.join("Chap5","sideband_postselected_phaseshift.json"),"5_2_phaseshift.pgf","Chap5_Fig5.2",size=1.0,cs=['b','k','r'],legend=True,h_pad=0,w_pad=0.,lloc=1)
#    TP.addPlot(os.path.join("Chap5","spectrum.json"),"5_1_spectrum.pgf","Chap5_Fig5.2",size=1.0,cs=['b','r','r','b'],legend=True,h_pad=0,w_pad=2.,lloc=1,tl=True)
#    TP.addPlot(os.path.join("Chap5","pol_spectra.json"),"5_4_spectra.pgf","Chap5_Fig5.4",size=1.0,cs=['b','b','r','k','r'],legend=False,tl=True,h_pad=0,w_pad=1.8,lloc=1,ls=['-','','-','-',''],m=['o','','o','',''])
#    TP.addPlot(os.path.join("Chap5","cond_phase_vs_density.json"),"5_7_phase.pgf","Chap5_Fig5.7",size=1.0,cs=['b','b','r','r'],legend=False,tl=False,h_pad=0,w_pad=1.8,lloc=1,ls=['-','','-',''],m=['o','','o','',''])
#    TP.addPlot(os.path.join("Chap5","propagation.json"),"5_5_propagation.pgf","Chap5_Fig5.5",size=1.0,cs=['b','b','r','r'],legend=False,tl=True,h_pad=0,w_pad=1.8,lloc=1,ls=['-','','-',''],m=['o','','o','',''])
    TP.addPlot(os.path.join("Chap5","storage_retrieval.json"),"5_6_storage.pgf","Chap5_Fig5.6",size=1.0,cs=['b','r'],legend=False,tl=False,h_pad=0,w_pad=0,lloc=1,ls=['','-'])
    
    TP.generatePlots()