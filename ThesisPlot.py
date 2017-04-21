import locale
import numpy as np
import io
import json
import pandas as pd
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
                        
                        ax.errorbar(x,y,yerr=yerr,label=self.dicts[d]['json'][sp_ax][sp]['label'])
                    elif self.dicts[d]['json'][sp_ax][sp]['type']=='plot':
                        df = pd.DataFrame.from_dict(json.loads(self.dicts[d]['json'][sp_ax][sp]['y']),orient='index')
                        df.set_index(np.array(df.index.values,dtype=np.float32),inplace=True)
                        df.sort(inplace=True)
                        x=np.array(df.index.values,dtype=np.float32)
                        y=np.array(df.values,dtype=np.float32)
                        
                        ax.plot(x,y,label=self.dicts[d]['json'][sp_ax][sp]['label'])
                    elif self.dicts[d]['json'][sp_ax][sp]['type']=='axh':
                        ax.axhline(np.float(self.dicts[d]['json'][sp_ax][sp]['y']),label=self.dicts[d]['json'][sp_ax][sp]['label'])
                    elif self.dicts[d]['json'][sp_ax][sp]['type']=='axv':
                        ax.axvline(np.float(self.dicts[d]['json'][sp_ax][sp]['y']),label=self.dicts[d]['json'][sp_ax][sp]['label'])
                        
                    ax.legend()
            
#            self.f.clear()
    
    def parsePlotDict(self,filename):
        with io.open(filename, 'r', encoding='utf-8') as f:
            plotDict=json.load(f)
            
        return plotDict

    def addPlot(self,name,outname,figid):
        self.dicts.update({figid:{'infile':name,'outfile':outname,'json':self.parsePlotDict(name)}})

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
    TP.addPlot("Chap2\Groupdelay\groupdelay.json","2_2_groupdelay.pgf","Chap2_Fig2.2")
#    TP.addPlot("Chap2\Groupdelay\groupdelay.json","2_3_groupdelay.pgf","Chap2_Fig2.3")
    
    TP.generatePlots()