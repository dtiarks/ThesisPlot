import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as c

a0=c.physical_constants["Bohr radius"][0]
kB=c.k
m=87*c.physical_constants["atomic mass constant"][0]
lD1 = 780*10**-9

A2nd=(3.28)

profile1=np.loadtxt("img_004/profileavg.had")

ref_avg=np.loadtxt("img_002/profileavg.had")


f=plt.figure(0)

#f.suptitle("Transversal integriertes Wolkenprofil")

ax0=f.add_subplot(111)

bckg_list= np.empty(220)
bckg_list.fill(455.0)
prof1=profile1[:,1]-ref_avg[:,1]
#prof1=profile1[:,1]-bckg_list


ax0.set_title("Longitudinal Density Profile", fontsize = 20)
ax0.plot(A2nd*profile1[:,0],prof1)

ax0.set_ylabel("OD (Arb. Un.)", fontsize = 18)
ax0.set_xlabel("Longitudinal Position ($\mu m$)", fontsize = 18)
#ax0.grid(True)
ax0.tick_params(axis='both', which='major', labelsize=16, width = 2, size = 5)

ax0.annotate("97$\mu m$", xy=(100, 150), xytext=(130, 145),
    arrowprops=dict(arrowstyle="->"), fontsize = 18)
ax0.annotate("", xy=(195, 150), xytext=(168, 150),
    arrowprops=dict(arrowstyle="->"))
    
plt.rcParams['axes.linewidth'] = 2

f.tight_layout()

f.show()
#f.savefig("PlugsProfile3.0MHz_Poster.png", transparent = True)