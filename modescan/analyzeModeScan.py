#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This script calls a variety of subroutines to analyze data from an OMC modescan
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
from get_PZT_sweeps import *
from fit_OMC_scan import *
#from fit_PRMI_scan import *
from calc_modescan import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import matplotlib.gridspec as gridspec
import pylab
matplotlib.rcParams.update({'savefig.dpi':150})

# functions called by this script:
#
# get_PZT_sweeps - scans data from OMC-PZT2_MON_DC_OUT and finds the times when the PZT was ramping
#                  also checks that peak heights make sense??
# 
# fit_OMC_scan   - fits modes to single PZT sweep, returns mode indices, heights, voltage to frequency fit, etc
#                  makes some assumptions about peak structure?  i.e. 45MHz sidebands
#
# plot_OMC_scan  - plots a single sweep, in terms of optical frequency
#                  plots peak fits and frequency fit
#
# calc_modescan  - calculates various interesting parameters using the modescan output
#                  mode matching for carrier & sidebands, alignment, contrast defect
#
# beacon_scan    - uses a known (large) excitation in DARM to measure coupling of arm carrier light to OMC

# ideally we are going to do this on the cluster, and grab lots of DCPD data at 16384Hz and PZT MON_DC data at 512Hz
# will need to interpolate PZT data

# These parameters need to be the same as used in fit_OMC_scan
FSR = 261.72e6
#HOMs = 57.39e6
HOMs = 0.21946*FSR

f1 = 9.10023e6
f2 = 5*f1

xmin = -60
xmax = 261+45+15


input_power = 2.8
a = input_power/2.8 # factor for pretty plotting

# for now we'll assume we have a ASCII output from DTT with matching sample rates

# string that says whether DARM offset was on
data_string = 'hipow'
#data = genfromtxt('data/full_lock_offset_' + data_string + '.txt')
#data = genfromtxt('data/full_lock_offset_' + data_string + '_27Aug2015.txt')
#data = genfromtxt('data/full_lock_offset_' + data_string + '_31Aug2015.txt')
#data = genfromtxt('data/PRMI_on_31Aug2015.txt')
data = genfromtxt('data/full_lock_' + data_string + '_5Sep2015_full.txt')

# note that times could be provided independently, using sample rate of time series
times = data[:,0]
pzt = data[:,1]
dcpd = data[:,2]

ymin = 1e-3
ymax = 2*max(dcpd)

fs = int(len(times)/(max(times) - min(times)))

# on the cluster we'll specify sweep times?
# for a txt file we'll specify indices

ramp_indices = get_PZT_sweeps(times,pzt)
sweeps = len(ramp_indices)

print
print 'We found', sweeps, 'modescans in this data.'

 
fignum=0

fignum=fignum+1
pylab.figure(fignum)

num_modes = 50
#num_modes = 32

# initialize some arrays to store calculations for each sweep
SB45_imbalance = zeros(sweeps)
CR_modematch = zeros(sweeps)
SB45_modematch = zeros(sweeps)
CR_align = zeros(sweeps)
SB45_align = zeros(sweeps)
CR_contrast = zeros(sweeps)
mode_fit_frequency = zeros((num_modes,sweeps))
mode_height = zeros((num_modes,sweeps))

for i in range(len(ramp_indices)):

    seg = ramp_indices[i]

    print
    print 'Fitting peaks in Segment', i
    print 'PZT start:', pzt[seg[0]], ', PZT stop:', pzt[seg[1]]

    pzt_scan = pzt[seg[0]:seg[1]]
    dcpd_scan = dcpd[seg[0]:seg[1]]

    mode_fits, mode_errors, volts_to_MHz, MHz_to_volts, mode_freqs, mode_text = fit_OMC_scan(pzt_scan, dcpd_scan, fs, fs, FSR, HOMs, f1, f2)

    mode_height[:,i] = mode_fits[:,0]

    peaks, params = shape(mode_fits)

    print 'We fit', peaks, 'modes.'

    # Problem: the modes are fit in terms of PZT voltage, but this is known to be nonlinear.  Lorentzian fits might be slightly assymetric.
    # It would be better to fit them in terms of optical frequency but there's a chicken-and-egg problem.

    #mode_fit_frequency[:,i], CR_contrast[i] = calc_modescan(mode_fits, mode_errors, volts_to_MHz, mode_freqs, mode_text, FSR)

    #for j in range(len(mode_text)):
    #    print mode_text[j], mode_fits[j,1], mode_errors[j,1], abs(mode_errors[j,1]*volts_to_MHz[1])

    optical_frequency = volts_to_MHz[0] + volts_to_MHz[1]*pzt_scan + volts_to_MHz[2]*pzt_scan**2 + volts_to_MHz[3]*pzt_scan**3
    pylab.semilogy(optical_frequency/1e6,dcpd_scan,'-',linewidth=0.8, alpha=0.8)



ymin = 2e-2
ymax = 2*max(dcpd_scan)

# The modes labeled here are just for illustration and do not reflect which modes were fit

#lsb_text = ['LSB1', 'LSB2', 'LSB1', 'LSB3', 'LSB4']
#lsb_locs = array([HOMs-f2, 2*HOMs-f2, HOMs+FSR-f2, 3*HOMs-f2, 4*HOMs-f2])

lsb_text = ['LSB1', 'LSB2', 'LSB1', 'LSB3', 'LSB4', 'LSB5', 'LSB5', 'LSB7']
lsb_locs = array([HOMs-f2, 2*HOMs-f2, HOMs+FSR-f2, 3*HOMs-f2, 4*HOMs-f2, 5*HOMs-f2, 5*HOMs-FSR-f2, 7*HOMs-FSR-f2])

usb_text = ['USB1', 'USB2', 'usb10', 'lsb9', 'lsb9']
usb_locs = array([HOMs+f2, 2*HOMs+f2, 10*HOMs-2*FSR+f1, 9*HOMs-FSR-f1, 9*HOMs-2*FSR-f1])

hom_text = ['CR1','CR2','CR3','CR4','CR5','CR6','CR7','CR8','CR9','CR10','CR11','CR12','CR13','CR14','CR15']
Hnum = 9

ehom_text = ['CR4','CR5','CR9','CR10','CR11','CR12','CR13']
ehom_locs = array([4*HOMs-FSR, 5*HOMs, 9*HOMs-2*FSR, 10*HOMs-2*FSR, 11*HOMs-2*FSR, 12*HOMs-2*FSR, 13*HOMs-2*FSR])

usb9_text = ['usb0','usb0','usb1', 'usb2','usb3','usb4','usb4','usb5','usb5','usb6','usb7','usb8']
usb9_locs = array([f1, FSR+f1,HOMs+f1, 2*HOMs+f1, 3*HOMs+f1, 4*HOMs+f1, 4*HOMs-FSR+f1, 5*HOMs+f1, 5*HOMs-FSR+f1, 6*HOMs-FSR+f1, 7*HOMs-FSR+f1, 8*HOMs-FSR+f1])

lsb9_text = ['lsb0', 'lsb1', 'lsb2', 'lsb3', 'lsb4', 'lsb4', 'lsb5', 'lsb5', 'lsb6', 'lsb7', 'lsb8', 'lsb0']
lsb9_locs = array([-f1, HOMs-f1, 2*HOMs-f1, 3*HOMs-f1, 4*HOMs-f1, 4*HOMs-FSR-f1, 5*HOMs-f1, 5*HOMs-FSR-f1, 6*HOMs-FSR-f1, 7*HOMs-FSR-f1, 8*HOMs-FSR-f1, FSR-f1])


# log spacing for plot labels
x = 10**arange(0.0,log10(a*0.75*ymax),log10(a*0.7*ymax)/6)
x = 10**arange(0.4,log10(0.75*ymax),log10(0.7*ymax)/12)

for i in range(2):
    pylab.plot([i*FSR/1e6, i*FSR/1e6],[ymin,ymax],'k-',linewidth=0.6)
    pylab.annotate('CR0', xy = (i*FSR/1e6, ymax), xytext = (1, -4), textcoords = 'offset points', ha = 'left', va = 'top',fontsize=10)

sb_text = ['USB0','LSB0','USB0','LSB0']
freq = array([FSR+f2, FSR-f2, f2, -1*f2])
for i in range(4):
    pylab.plot([freq[i]/1e6, freq[i]/1e6],[ymin,ymax],'r-',linewidth=0.6)
    pylab.annotate(sb_text[i], xy = (freq[i]/1e6, ymax), xytext = (1, -10*a), textcoords = 'offset points', ha = 'left', va = 'top',fontsize=10)

for i in range(len(usb9_locs)):
    pylab.plot([usb9_locs[i]/1e6, usb9_locs[i]/1e6],[ymin,x[5]*a],'m-',linewidth=0.4)
    pylab.annotate(usb9_text[i], xy = (usb9_locs[i]/1e6, x[5]*a), xytext = (1, 0), textcoords = 'offset points', ha = 'left', va = 'bottom',fontsize=7)

for i in range(len(lsb9_locs)):
    pylab.plot([lsb9_locs[i]/1e6, lsb9_locs[i]/1e6],[ymin,x[4]*a],'m-',linewidth=0.4)
    pylab.annotate(lsb9_text[i], xy = (lsb9_locs[i]/1e6, x[4]*a), xytext = (1, 0), textcoords = 'offset points', ha = 'left', va = 'bottom',fontsize=7)

for i in range(len(lsb_locs)):
    pylab.plot([lsb_locs[i]/1e6, lsb_locs[i]/1e6],[ymin,x[3]*a],'g-',linewidth=0.4)
    pylab.annotate(lsb_text[i], xy = (lsb_locs[i]/1e6, x[3]*a), xytext = (1, 0), textcoords = 'offset points', ha = 'left', va = 'bottom',fontsize=7)

for i in range(len(usb_locs)):
    pylab.plot([usb_locs[i]/1e6, usb_locs[i]/1e6],[ymin,x[2]*a],'g-',linewidth=0.4)
    pylab.annotate(usb_text[i], xy = (usb_locs[i]/1e6, x[2]*a), xytext = (1, 0), textcoords = 'offset points', ha = 'left', va = 'bottom',fontsize=7)

for i in range(Hnum):
    pylab.plot([remainder((i+1)*HOMs,FSR)/1e6, remainder((i+1)*HOMs,FSR)/1e6],[ymin,x[1]*a],'b-',linewidth=0.4)
    pylab.annotate(hom_text[i], xy = (remainder((i+1)*HOMs,FSR)/1e6, x[1]*a), xytext = (1, 0), textcoords = 'offset points', ha = 'left', va = 'bottom',fontsize=7)

for i in range(len(ehom_locs)):
    pylab.plot([ehom_locs[i]/1e6, ehom_locs[i]/1e6],[ymin,x[0]*a],'b-',linewidth=0.4)
    pylab.annotate(ehom_text[i], xy = (ehom_locs[i]/1e6, x[0]*a), xytext = (1, 0), textcoords = 'offset points', ha = 'left', va = 'bottom',fontsize=7)

pylab.grid(True, which='both', linestyle=':', alpha=0.4)

pylab.ylabel('OMC DCPD Sum [mA]')
pylab.xlabel('Optical Frequency [MHz]')

pylab.xlim(xmin,xmax)
pylab.ylim(ymin,ymax)

pylab.savefig('OpFreq_fit_offset_' + data_string + '.pdf')





fignum=fignum+1
pylab.figure(fignum)

optical_frequency = volts_to_MHz[0] + volts_to_MHz[1]*pzt_scan + volts_to_MHz[2]*pzt_scan**2 + volts_to_MHz[3]*pzt_scan**3

pfit = zeros(shape(pzt_scan))
for i in range(peaks):
    pk = mode_fits[i,0]/(1.0+((pzt_scan-mode_fits[i,1])/mode_fits[i,2])**2)
    pfit += pk
    pylab.semilogy(optical_frequency/1e6,pk,'b-',linewidth=0.8)

pylab.semilogy(optical_frequency/1e6,dcpd_scan,'k-',linewidth=0.8, alpha=0.8)
pylab.semilogy(optical_frequency/1e6,pfit,'orange',linestyle='-',linewidth = 1.8,alpha=0.8)

pylab.grid(True, which='both', linestyle=':', alpha=0.4)

pylab.ylabel('OMC DCPD Sum [mA]')
pylab.xlabel('Optical Frequency [MHz]')

pylab.xlim(xmin,xmax)
pylab.ylim(ymin,ymax)

pylab.savefig('test_fit_offset_' + data_string + '.pdf')



# Plot errors in frequency fit from the last sweep
# Should set this up to plot an arbitrary sweep
    

fignum=fignum+1
pylab.figure(fignum)

gs = gridspec.GridSpec(2,1,height_ratios=[3,1])

fit_frequencies = volts_to_MHz[0] + volts_to_MHz[1]*mode_fits[:,1] + volts_to_MHz[2]*mode_fits[:,1]**2 + volts_to_MHz[3]*mode_fits[:,1]**3

pylab.subplot(gs[0])

pylab.plot(pzt_scan,optical_frequency/1e6,'k-',linewidth=0.8)
pylab.plot(mode_fits[:,1],fit_frequencies/1e6,'bo')

pylab.grid(True, which='both', linestyle=':', alpha=0.4)
pylab.ylabel('Optical Frequency [MHz]')
#pylab.xlim(4,63)
pylab.ylim(xmin,xmax) # because x is usually frequency coordinate

pylab.subplot(gs[1])

pylab.errorbar(mode_fits[:,1],(fit_frequencies-mode_freqs)/1e6,yerr=std(mode_fit_frequency,1)/1e6,fmt='bo',ecolor='k',markersize=5)

"""
for i in range(len(mode_text)):
#    print mode_text[i], (fit_frequencies[i]-mode_freqs[i])/1e6
#    print mode_text[i], mode_freqs[i]/1e6, mode_fit_frequency[i,:]/1e6, std(mode_fit_frequency[i,:])/1e6
    print mode_text[i], mode_fits[i,0], mode_freqs[i]/1e6, std(mode_fit_frequency[i,:])/1e6
#    print mode_text[i], mode_fits[i,0]
"""

pylab.grid(True, which='both', linestyle=':', alpha=0.4)
pylab.xlabel('PZT Output [V]')
pylab.ylabel('residuals [MHz]')
#pylab.xlim(4,63)

pylab.savefig('freq_fit_offset_' + data_string + '.pdf')



# Plot mode height across many sweeps

fignum=fignum+1
pylab.figure(fignum)

optical_frequency = volts_to_MHz[0] + volts_to_MHz[1]*pzt_scan + volts_to_MHz[2]*pzt_scan**2 + volts_to_MHz[3]*pzt_scan**3

pfit = zeros(shape(pzt_scan))
for i in range(peaks):
    pk = mode_fits[i,0]/(1.0+((pzt_scan-mode_fits[i,1])/mode_fits[i,2])**2)
    pfit += pk
    pylab.semilogy(optical_frequency/1e6,pk,'b-',linewidth=0.8)

pylab.semilogy(optical_frequency/1e6,dcpd_scan,'k-',linewidth=0.8, alpha=0.8)
pylab.semilogy(optical_frequency/1e6,pfit,'orange',linestyle='-',linewidth = 1.8,alpha=0.8)

pylab.grid(True, which='both', linestyle=':', alpha=0.4)

pylab.ylabel('OMC DCPD Sum [mA]')
pylab.xlabel('Optical Frequency [MHz]')

pylab.xlim(xmin,xmax)
pylab.ylim(ymin,ymax)

pylab.savefig('test_fit_offset_' + data_string + '.pdf')





# The following code is for plotting the frames of the movie

pidx = arange(0,len(mode_text),1)


dcpd_scan_new = dcpd_scan
for i in range(peaks):

    fignum=fignum+1
    pylab.figure(fignum)

    pylab.semilogy(optical_frequency/1e6,dcpd_scan_new,'k-',linewidth=0.8)

    pk = mode_fits[pidx[i],0]/(1.0+((pzt_scan-mode_fits[pidx[i],1])/mode_fits[pidx[i],2])**2)
    pylab.semilogy(optical_frequency/1e6,pk,'orange',linestyle='-',linewidth = 1.8,alpha=0.8)

    pylab.annotate(mode_text[pidx[i]], xy = (120, 30), xytext = (1, 0), textcoords = 'offset points', ha = 'left', va = 'bottom',fontsize=16)

    dcpd_scan_new = dcpd_scan_new - pk

    pylab.grid(True, which='both', linestyle=':', alpha=0.4)

    pylab.ylabel('OMC DCPD Sum [mA]')
    pylab.xlabel('Optical Frequency [MHz]')
    pylab.xlim(xmin,xmax)
    pylab.ylim(ymin,ymax)

    if i < 10:
        pylab.savefig('images/peak_fit_0' + str(i) + '.png')
    else:
        pylab.savefig('images/peak_fit_' + str(i) + '.png')

    pylab.close()


