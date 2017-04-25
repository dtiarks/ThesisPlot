import locale
import numpy as np
import io
import json
import pandas as pd
import ast
# Set to German locale to get comma decimal separater
locale.setlocale(locale.LC_NUMERIC, 'deu_deu')
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
    "figure.figsize": [5.67, 1.75],   # default fig size of 0.9 textwidth
    "errorbar.capsize": 0,             # set standard
    "markers.fillstyle": 'none',
    "lines.markersize": 4,
    "lines.linewidth": 1.5,
    "legend.fancybox": True,
    "text.latex.preamble": preamble,
    #    "pgf.debug" : True,
    "pgf.preamble": preamble,
    "legend.numpoints": 1,
    "axes.formatter.use_locale": True
}

mpl.rcParams.update(pgf_with_latex)

import matplotlib.pyplot as plt


class ThesisPlot(object):

    def __init__(self):
        self.dicts=dict()

    def generatePlots(self):
        for d in self.dicts:
            self.f=plt.figure()
            
            for sp_ax in self.dicts[d]['json']:
                ax=self.f.add_subplot(sp_ax)
                
                for sp in self.dicts[d]['json'][sp_ax]:
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
                        
                        ax.errorbar(x,y,yerr=yerr,label=label)
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
                        
                        ax.plot(x,y,label=label)
                    elif self.dicts[d]['json'][sp_ax][sp]['type']=='axh':
                        ax.axhline(np.float(self.dicts[d]['json'][sp_ax][sp]['y']),label=self.dicts[d]['json'][sp_ax][sp]['label'])
                    elif self.dicts[d]['json'][sp_ax][sp]['type']=='axv':
                        ax.axvline(np.float(self.dicts[d]['json'][sp_ax][sp]['y']),label=self.dicts[d]['json'][sp_ax][sp]['label'])
                        
#                    ax.legend()
                    
            if len(self.f.axes) > 1:
                for n, ax in enumerate(self.f.axes):
                    ax.text(0.1, 0.9, r'\textbf{(' + chr(len(self.f.axes)-1-n + 97) + ')}', transform=ax.transAxes,
                            weight='bold', ha='center', va='center')
                    # label position in window coordinates
                    xwin, ywin = ax.transAxes.transform((0.1, 0.9))
                    for l in ax.yaxis.get_major_ticks():
                        # check if a label overlaps with enumeration
                        bbox = l.label1.get_window_extent()
                        print bbox, xwin, ywin
                        if self._overlaps(np.array(bbox), xwin, ywin):
                            l.label1.set_visible(False)
            
            s=self.figsize(self.dicts[d]['size'],1.0)
            self.f.set_size_inches(*s)
            print self.dicts[d]['outfile']
            self.f.savefig(self.dicts[d]['outfile'])
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

    def addPlot(self,name,outname,figid,size=2):
        self.dicts.update({figid:{'infile':name,'outfile':outname,'size':size,'json':self.parsePlotDict(name)}})

    def figsize(self, rows, scale):
        # Get this from LaTeX using \the\textwidth
        fig_width_pt = 455.24416
        inches_per_pt = 1.0 / 72.27                       # Convert pt to inch
        # Aesthetic ratio (you could change this) * 0.5
        golden_mean = rows * 0.5 * (np.sqrt(5.0) - 1.0) / 2.0
        fig_width = fig_width_pt * inches_per_pt * scale    # width in inches
        fig_height = fig_width * golden_mean              # height in inches
        fig_size = [fig_width, fig_height]
        return fig_size

if __name__=='__main__':
    TP=ThesisPlot()
#    TP.addPlot("Chap2\Groupdelay\groupdelay.json","2_2_groupdelay.pgf","Chap2_Fig2.2")
    TP.addPlot("Chap2\Suszept\suszept.json","2_1_suszept.pgf","Chap2_Fig2.1",size=2)
    
    TP.generatePlots()