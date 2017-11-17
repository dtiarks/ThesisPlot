# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 17:49:13 2016

@author: tstolz
"""

import numpy as np
from scipy.special import erf
PHI = lambda x: 0.5 * (1 + erf(x / np.sqrt(2)))
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.animation as animation
plt.rcParams['animation.ffmpeg_path'] = 'ffmpeg'
plt.ion()
import time
import json
import pandas as pd
import io

### PARAMETERS AND CONSTANTS

hbar = 6.626070040e-34/(2 * np.pi) # Js, Planck constant, CODATA 2014
d = 2.534e-29 # Cm, dipole matrix element (D. A. Steck) of the 87Rb D2 Cycling Transition
Gamma_e = 2*np.pi * 5.746e6 # decay rate (D. A. Steck)
epsilon_0 = 8.854187817e-12 # dielectric constant, CODATA 2014
L = 48e-6 # medium length in m
omega_s = 2*np.pi * 377.23e12 # rad/s, transition frequency
c = 299792458 # m/s, speed of light CODATA 2014

# SIMULATION INPUT

# simulation parameters
sim_pars = {
                "tL" : 300e-6, # length of the time interval to include. determines the frequency resolution.
                "N" : 2**18, # number of points to include. determines the time resolution. powers of two yield best FFT efficiency

                # experiment parameters
                "rho_peak" : 0.5*3.1e12/1e-6, # peak density in m^-3 (cm^-3/centi^-3)
                "L" : 48e-6, # medium length in m
          
                "gs_EIT" : {
                            "Delta_c" : 0, #2*np.pi * 100e3, # detuning of the coupling laser
                            "gamma_21" : 1, #2*np.pi*300e3, # dephasing rate (rad/s)
                            "Omega_c" : 0, #2*np.pi*2e6, #3e6 # coupling Rabi-Frequency (rad/s)
                            "ladder" : True,
                        },   
                "ryd_EIT" : {
                            #"Delta_c" : -1.81*Gamma_e, #-50e3, # detuning of the coupling laser
                            "Delta_c" : -1.4*Gamma_e, #-50e3, # detuning of the coupling laser
                            "gamma_21" : 2*np.pi* 1100e3, #500e3, # dephasing rate (rad/s)
                            "Omega_c" : 0.8*Gamma_e,#10e6, # coupling Rabi-Frequency (rad/s)
                            "ladder" : False,
                         },
                         
                # input pulse
                "tRMS" : 0.01e-6, # RMS width of the input gaussian
                "tW" : 0.6e-6, # width of the rectangular input pulse
                "Delta_s0" : -1.4*Gamma_e, # detuning of the signal laser from resonance (rad/s)
                "S3_in" : -0.75, # Input S3 parameter, translates to (p_R - p_L) / (p_R + p_L)
                "in_func" :  '''lambda times, **kwargs: np.sqrt(  PHI((kwargs["tL"]/2 - kwargs["tW"]/2 - times)/(-0.0005e-6) ) * PHI((kwargs["tL"]/2 + kwargs["tW"]/2 - times)/0.0005e-6) ) * np.exp(-1j * kwargs["Delta_s0"] * times)'''#'''lambda times, **kwargs: np.sqrt( np.exp(- 0.5 * (times - kwargs["tL"]/2)**2/(kwargs["tRMS"])**2) * PHI((kwargs["tL"]/2 - times)/0.015e-6) ) * np.exp(-1j * kwargs["Delta_s0"] * times)'''
                # the half-gaussian envelope of the pulse as in the experiment control code. using an error function (PHI) to smoothen
                # the cut like the AOM would do (with a time constant of about 20 ns). 
            }

SHOW_RESULT = False
SAVE_ANIMATION = False



# utility functions

def susceptibility(Delta_s, Delta_c, gamma_21, Omega_c, ladder=True):
    delta = (Delta_s + (-1 + 2*int(ladder)) * Delta_c) # two photon detuning
    return 1j*(gamma_21 - 2j * delta)/(np.abs(Omega_c)**2 + (Gamma_e - 2j * Delta_s)*(gamma_21 - 2j * delta))

def norm(envelope, tL, N):
    Dt = tL/float(N)
    return np.sum(np.abs(envelope)**2)*Dt

def getSpectrum(envelope, tL, N):
    spectrum = tL * np.fft.ifft(envelope)
    ordered_spectrum = spectrum[range(N/2, N) + range(N/2)] # the negative frequencies are in the second half of the spectrum
    deltas = 2*np.pi*(np.arange(N)-N/2)/tL
    return deltas, ordered_spectrum

def propagateSpectrum(deltas, spectrum, L, chi):
    #return spectrum*np.exp(1j*omega_s/c * L * np.sqrt(1 + chi(deltas)))
    return spectrum*np.exp(1j*omega_s/c * L * (1 + 0.5*chi(deltas)))
    
def getEnvelope(spectrum, tL, N):
    reordered_spectrum = spectrum[range(N/2, N) + range(N/2)] # reset to standard order
    envelope = np.fft.fft(reordered_spectrum)/tL
    times = np.arange(N)*tL/N
    return times, envelope
       
# the propagation routine
def simulate_propagation(sim_pars):
    start_time = time.time()
#    f = open("EIT_Propagation_Parameters_%d.json" % start_time, "wb")
#    json.dump( sim_pars, f, indent=4)
#    f.close()
    K = sim_pars
    tL, N, rho_peak, L, S3_in, gs_EIT, ryd_EIT = [ K["tL"], K["N"], K["rho_peak"], K["L"], K["S3_in"], K["gs_EIT"], K["ryd_EIT"]]
    times = np.arange(N)*tL/N # array containing the sampling times
    
    
    env_non_norm = eval(K["in_func"])(times, **K) # in_func is a lambda function string
    envelope = env_non_norm/np.sqrt(norm(env_non_norm, tL, N))
    
    
    deltas, spectrum = getSpectrum(envelope, tL, N)
    
    chi_0 = 2*rho_peak*d**2 / (epsilon_0*hbar*Gamma_e) # prefactor of the susceptibility for the Cycling transition |2,2> -> |2',3'>
    chi_gs = lambda d:  0*chi_0*Gamma_e*susceptibility(d, **gs_EIT)/6 # transition strength is 1/6
    chi_ryd = lambda d:  1*chi_0*Gamma_e*susceptibility(d, **ryd_EIT) 
    
    print "%f"%(chi_0)
    
    spec_at_L_gs = propagateSpectrum(deltas, spectrum, L, chi_gs) 
    spec_at_L_ryd = propagateSpectrum(deltas, spectrum, L, chi_ryd)
    
    times_gs, envelope_gs = getEnvelope(spec_at_L_gs, tL, N)
    times_ryd, envelope_ryd = getEnvelope(spec_at_L_ryd, tL, N)
    phases = (np.angle(envelope_gs)-np.angle(envelope_ryd)) % (2*np.pi) - np.pi
    
    int_gs = 0.5*(1+S3_in)*np.abs(envelope_gs)**2
    int_ryd = 0.5*(1-S3_in)*np.abs(envelope_ryd)**2
    S3 = (int_gs - int_ryd)/(int_gs + int_ryd)
#    S3 = (envelope_gs - envelope_ryd) / (envelope_gs + envelope_ryd)
    S2 = np.sin(phases) * np.sqrt( 1 - S3**2 )
    S1 = np.cos(phases) * np.sqrt( 1 - S3**2 )
    print "finished after %f seconds" % (time.time()-start_time)
    
    microsecs = 1e6*(times - tL/2)
    xaxis=(Gamma_e*(times - tL/2)+57)
    print xaxis[len(microsecs)/2-2400:len(microsecs)/2+2400]
    
    pulse_in=np.abs(envelope)**2
    max_in=np.max(pulse_in)
    pulse_in/=max_in
    pulse_out=np.abs(envelope_ryd)**2
    pulse_out/=max_in
    
    plot_dict={}
    
    plt.figure(0)
    
    plt.subplot(121)
#    plt.plot(microsecs, 0.5*(1+S3_in)*np.abs(envelope_gs)**2, 'k-', label='GS, transm. = %2.3f %%' % (100*norm(envelope_gs, tL, N)), color='orange')
    plt.plot(xaxis[len(microsecs)/2-2400:len(microsecs)/2+2400], pulse_out[len(microsecs)/2-2400:len(microsecs)/2+2400], 'k-', label='Ryd, transm. = %2.3f %%' % (100*norm(envelope_ryd, tL, N)), color='cyan')
    plt.plot(xaxis[len(microsecs)/2-2400:len(microsecs)/2+2400], pulse_in[len(microsecs)/2-2400:len(microsecs)/2+2400])
    plt.xlim(-57.,171.)
    
    h=pd.DataFrame(index=xaxis[len(microsecs)/2-2400:len(microsecs)/2+2400],data=pulse_out[len(microsecs)/2-2400:len(microsecs)/2+2400])
    h2=pd.DataFrame(index=xaxis[len(microsecs)/2-2400:len(microsecs)/2+2400],data=pulse_in[len(microsecs)/2-2400:len(microsecs)/2+2400])
    plot_dict['121']={
        'A':{'type':'plot','y':h2[0].to_json(),'xlim':(-57.,171.)},
        'B':{'type':'plot','y':h[0].to_json(),'ylabel':u'Normierte Intensit\"at\quad','xlabel':u'Zeit ($1/\Gamma_3$)','xlim':(-57.,171.),'ylim':(0.,1.1),'num':'a'}
        
    }
    
    plt.subplot(122)
    nlow=1250
    nhigh=1400
    nhigh=nlow
    plt.plot(xaxis[len(microsecs)/2-nlow:len(microsecs)/2+nhigh], phases[len(microsecs)/2-nlow:len(microsecs)/2+nhigh], color='purple')
#    plt.axhline(np.mean(phases_cw[len(microsecs)/2-1400:len(microsecs)/2+1400]))
    plt.axhline(-0.5*omega_s/c * L *chi_0*Gamma_e*np.real(susceptibility(K["Delta_s0"], **ryd_EIT) )%(2*np.pi)-np.pi)
    plt.xlim(-57.,171.)
    plt.ylim(-np.pi, np.pi)
    
    h=pd.DataFrame(index=xaxis[len(microsecs)/2-nlow:len(microsecs)/2+nhigh],data=phases[len(microsecs)/2-nlow:len(microsecs)/2+nhigh])
    plot_dict['122']={
        'A':{'type':'axh','y':-0.5*omega_s/c * L *chi_0*Gamma_e*np.real(susceptibility(K["Delta_s0"], **ryd_EIT) )%(2*np.pi)-np.pi},
        'B':{'type':'plot','y':h[0].to_json(),'ylabel':u'Phase (rad)','xlabel':u'Zeit ($1/\Gamma_3$)','xlim':(-57.,171.),'ylim':(-1*np.pi,np.pi),'num':'b'}
        
    }
    
#    with io.open('eit_propagation.json', 'w+') as f:
#        f.write(json.dumps(plot_dict, ensure_ascii=False,indent=4))
       
    if SHOW_RESULT:
        
        CHI_gs = chi_gs(deltas)
        CHI_ryd = chi_ryd(deltas)
        
        plt.figure("EIT Propagation", figsize=(18,10))
        
        plt.subplot(331)
        plt.title('Input Signal')
        plt.plot(microsecs, np.abs(envelope)**2, 'k-')
        plt.xlabel("Time ($\mu$s)")
        plt.ylabel("Intensity (a.u.)")
        plt.xlim( - 3, + 3)
        plt.grid()
        
        plt.subplot(334)
        plt.plot(microsecs, np.angle(envelope), 'r-')
        plt.xlabel("Time ($\mu$s)")
        plt.ylabel("Phase (rad)")
        plt.grid()
        plt.xlim(- 3, + 3)
        
        plt.subplot(337)
        plt.plot(deltas/(2e6*np.pi), np.abs(spectrum)**2, 'k-')
        plt.xlabel("Frequency (MHz)")
        plt.ylabel("Power Spectrum (a.u.)")
        plt.xlim(-20,20)
        plt.grid()
        
        plt.subplot(332)
        plt.title("Propagating with Susceptibility $\chi$")
        plt.plot(deltas/(2e6*np.pi), np.imag(CHI_gs), color='orange', label=r'Re($\chi_{GS}$)')
        plt.plot(deltas/(2e6*np.pi), np.real(CHI_gs), color='red', label='Im($\chi_{GS}$)')
        plt.plot(deltas/(2e6*np.pi), np.imag(CHI_ryd), color='cyan', label=r'Re($\chi_{Ryd}$)')
        plt.plot(deltas/(2e6*np.pi), np.real(CHI_ryd), color='blue', label='Im($\chi_{Ryd}$)')
        plt.xlabel("Frequency (MHz)")
        plt.xlim(-20,20)
        leg = plt.legend(fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.grid()
        
        plt.subplot(335)
        plt.plot(deltas/(2e6*np.pi), np.exp(- omega_s/c * L * np.imag(CHI_gs)), color='orange', label='GS')
        plt.plot(deltas/(2e6*np.pi), np.exp(- omega_s/c * L * np.imag(CHI_ryd)), color='cyan', label='Ryd')
        plt.ylabel('Transmission')
        plt.xlabel("Frequency (MHz)")
        plt.xlim(-20,20)
        plt.grid()
        plt.axvline(K["Delta_s0"]/(2e6*np.pi))
        leg = plt.legend(fancybox=True)
        leg.get_frame().set_alpha(0.5)
        
        plt.subplot(338)
        plt.plot(deltas/(2e6*np.pi), omega_s/(2*c) * L *np.real(CHI_gs), color='red', label='GS')
        plt.plot(deltas/(2e6*np.pi), omega_s/(2*c) * L *np.real(CHI_ryd), color='blue', label='Ryd')
        leg = plt.legend(fancybox=True)
        leg.get_frame().set_alpha(0.5)
        plt.ylabel('Phase (rad)')
        plt.xlabel("Frequency (MHz)")
        plt.xlim(-20,20)
        plt.grid()
        
        plt.subplot(333)
        plt.title('Output Signal')
        plt.plot(microsecs, 0.5*(1+S3_in)*np.abs(envelope_gs)**2, 'k-', label='GS, transm. = %2.3f %%' % (100*norm(envelope_gs, tL, N)), color='orange')
        plt.plot(microsecs, 0.5*(1-S3_in)*np.abs(envelope_ryd)**2, 'k-', label='Ryd, transm. = %2.3f %%' % (100*norm(envelope_ryd, tL, N)), color='cyan')
        plt.xlabel("Time ($\mu$s)")
        plt.ylabel("Intensity (a.u.)")
        plt.xlim(- 3, 3)
        plt.grid()
        leg = plt.legend(fancybox=True)
        leg.get_frame().set_alpha(0.5)
        
        plt.subplot(336)
        plt.plot(microsecs, phases, color='purple')
        plt.xlabel("Time ($\mu$s)")
        plt.ylabel("Phase Difference (rad)")
        plt.grid()
        plt.xlim(- 3, + 3)
        plt.ylim(-np.pi, np.pi)
        
        plt.subplot(339)
        plt.plot(microsecs, S1, label="S1")
        plt.plot(microsecs, S2, label="S2")
        plt.plot(microsecs, S3, label="S3")
        plt.xlabel("Time ($\mu$s)")
        plt.grid()
        plt.xlim(- 3, + 3)
        plt.ylim(-1,1)
        leg = plt.legend(fancybox=True)
        leg.get_frame().set_alpha(0.5)
        
#        plt.subplot(339)
#        plt.plot(deltas/(2e6*np.pi), 0.5*(1+S3_in)*np.abs(spec_at_L_gs)**2, color='red', label='GS')
#        plt.plot(deltas/(2e6*np.pi), 0.5*(1-S3_in)*np.abs(spec_at_L_ryd)**2, color='blue', label='Ryd')
#        plt.xlabel("Frequency (MHz)")
#        plt.ylabel("Power Spectrum (a.u.)")
#        plt.grid()
#        plt.xlim(-5,5)
        
        plt.tight_layout()
        
        plt.savefig('results/EIT_propagation_result_%d.png' % start_time)
    return envelope_gs, envelope_ryd
    
def animate_propagation(slowdown=1e7, timespan = [-1,2], points=200, zrange=[-150,250], **kwargs):
    AMPS = []
    microsecs = 1e6*(times - tL/2)
    indices = np.where((timespan[0] < microsecs) & (microsecs < timespan[1]))[0]
    zvals = np.linspace(zrange[0], zrange[1], points)
    envelope = env/np.sqrt(norm(env, tL, N))
    
    for z in zvals[zvals<=0]: # free propagation; create a time shift by an index shift
        shift = -int((z*1e-6/c)/(tL/N)) # number of time intervals within the travelled time, negative because we have to go to earlier times
        AMPS.append(np.abs(envelope[indices+shift])**2)
        
    deltas, spectrum = getSpectrum(envelope, tL, N) # fourier transform for propagation (save time, only do it once)
    chi = lambda d:  chi_0*Gamma_e*susceptibility(d, **gs_EIT)
    
    for z in zvals[(zvals>0) & (zvals*1e-6<L)]: # propagate through medium until L
        spec_at_z = propagateSpectrum(deltas, spectrum, z*1e-6, chi)
        times_f, envelope_f = getEnvelope(spec_at_z, tL, N)
        AMPS.append(np.abs(envelope_f[indices])**2)
    spec_at_L = propagateSpectrum(deltas, spectrum, L, chi)
    times_f, envelope_f = getEnvelope(spec_at_L, tL, N) 
    for z in zvals[zvals*1e-6>=L]: # free propagation again
        shift = -int((z*1e-6/c)/(tL/N))
        AMPS.append(np.abs(envelope_f[indices+shift])**2) # need to take the last pulse shape now
    AMPS = np.transpose(AMPS) # get an array of position waveforms for each time instead of an array of time waveforms for each position

    def update_line(num, line): # function for the animation
        line.set_data(zvals, AMPS[num])
        return line,

    fig = plt.figure()
    l, = plt.plot(zvals, AMPS[0], 'r-')
    plt.xlim(*zrange)
    plt.ylim(0, np.max(AMPS))
    plt.xlabel('Position ($\mu$m)')
    currentAxis = plt.gca()
    currentAxis.add_patch(Rectangle((0,0), 60, np.max(AMPS), facecolor="grey"))
    ani = animation.FuncAnimation(fig, update_line, len(AMPS), fargs=[l], interval=tL/N*slowdown*1e3, blit=True)
    if SAVE_ANIMATION:
        ani.save('EIT_Propagation_%d.mp4' % time.time(), metadata={'artist':'tstolz'}, writer="ffmpeg")
    else:
        plt.show()

    
    
if __name__=='__main__':
    simulate_propagation(sim_pars)
    
#    slowdown = 1e7
#    timespan = [-1,2]
#    points=300
#    zrange= [-150, 250] #[-300e6,300e6]# [-150,250] #
#    
#    SAVE_ANIMATION=True
#    animate_propagation(slowdown, timespan, points, zrange)
