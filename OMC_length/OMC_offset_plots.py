#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import cmath as cm
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab

fig_width_pt = 600  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio, = 0.618
golden_mean = 0.6
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width,fig_height]
matplotlib.rcParams.update({'figure.figsize': fig_size})

lam = 1064.0e-9

# OMC input & output couplers are 8000ppm transmission
r_i = sqrt(1-8000e-6)
r_e = sqrt(1-8000e-6)

t_i = sqrt(1-r_i**2)
t_e = sqrt(1-r_e**2)

L = 0
finesse = pi*sqrt(r_i * r_e) / (1- r_i*r_e)



data1 = genfromtxt('data/OMC_meters_13Jun2015.txt')
freq1 = data1[:,0]
deltaL = data1[:,1]

data2 = genfromtxt('data/OMC_RIN_13Jun2015.txt')
freq2 = data2[:,0]
RIN = data2[:,1]
NULL = data2[:,2]


data3 = genfromtxt('data/OMC_offset_RIN_coh_13Jun2015.txt')
freq3 = data3[:,0]
ISS = data3[:,1]
POP = data3[:,2]
ASAIR = data3[:,3]

data4 = genfromtxt('data/OMC_offset_RIN_suppressed_13Jun2015.txt')
freq4 = data4[:,0]
RIN_low = data4[:,1]

data5 = genfromtxt('data/OMC_meters_suppressed_13Jun2015.txt')
freq5 = data5[:,0]
deltaL_low = data5[:,1]

data6 = genfromtxt('data/L1_data.txt',delimiter=",")
freq6 = data6[:,0]
deltaL_L1 = data6[:,1]

data7 = genfromtxt('data/H1_PDH_data.txt',delimiter=",")
freq7 = data7[:,0]
deltaL_H1_PDH = data7[:,1]

data8 = genfromtxt('data/OMC_DCPD_dark_RIN_13Jun2015.txt')
freq8 = data8[:,0]
dark = data8[:,1]

data9 = genfromtxt('data/CAL_DeltaL_7Jun.txt')
freq9 = data9[:,0]
CAL = data9[:,1]
DARM = data9[:,2]

data10 = genfromtxt('data/OMC_cavity_noise_projected.txt')
freq10 = data10[:,0]
OMC_x2_ASD = data10[:,1]


"""
# Calculate fluctuation in DCPD Sum as a function of frequency
# Use the suppressed length noise
dx = 1e7*deltaL_low
E_trans = ( t_i * t_e * exp(1.0j*2*pi*dx/lam))/(1 - r_i * r_e * exp(2.0j*2*pi*dx/lam))
P_trans = real(E_trans * conj(E_trans))
"""

# Calculate delta_P using an approximation for the OMC cavity resonance
f_rms = freq5
dx = deltaL_low

#f = freq7
#dx = deltaL_H1_PDH

#P_trans = 2.1e18 * (dx**2)

# Use RIN = (4*finesse/lambda * x)**2
# (4*finesse/lambda) = 2.16e18 for OMC, finesse=391
# This fits the formula for transmitted E-field very well for RMS motion of ~3e-11 meters - should be ok for suppressed motion
f = freq10
P_trans = 2.16e18*OMC_x2_ASD

# Use RIN = (4*finesse/lambda)**2 * x_rms * x
x_rms = 3e-13
P_trans_rms = 2.16e18 * 2 * x_rms * dx

# Convert P_trans to DARM_IN1 using OMC-READOUT calculation
Pref = 1.56
P0 = 24.1
xf = 14.0
x0 = 15.8
G = 8.555e-7
OMC_ERR = ((P_trans * Pref / P0)*(xf**2/2/x0)) * G
8
OMC_ERR_rms = ((P_trans_rms * Pref / P0)*(xf**2/2/x0)) * G

#OMC_DARM = OMC_ERR * 14e-12 / 1e-5
#OMC_DARM = OMC_ERR * 5e-7

# Very approximate model of DARM response - zero at 389Hz, pole at 7000Hz
# The constant factor calibrates DARM_IN1 counts --> meters.  This is 3x smaller than expected from DARM offset calibration...
CAL_TF = 5.0e-7 * abs((1 + 1.0j*f/389.0)/(1 + 1.0j*f/7000.0))

CAL_TF_rms = 5.0e-7 * abs((1 + 1.0j*f_rms/389.0)/(1 + 1.0j*f_rms/7000.0))

OMC_DARM = OMC_ERR * CAL_TF
OMC_DARM_rms = OMC_ERR_rms * CAL_TF_rms

print deltaL_low
print P_trans

matplotlib.rcParams.update({'savefig.dpi':250})

fignum=0


fignum=fignum+1
pylab.figure(fignum)

pylab.loglog(freq4,RIN_low,'r-',linewidth=0.8,label='OMC DCPD RIN, suppressed')
pylab.loglog(freq2,RIN,'b-',linewidth=0.8,label='OMC DCPD RIN, free-running')
pylab.loglog(freq8,dark,'-',color='0.6',linewidth=0.8,label=r'OMC DCPD Dark ($\div$ by SUM)')
pylab.loglog(freq2,NULL,'k-',linewidth=0.8,label=r'OMC DCPD NULL ($\div$ by SUM)')
pylab.loglog([min(freq2),max(freq2)],[5.6e-9, 5.6e-9],'y--',linewidth=1.4,label='Shot Noise, 10.2mA')
pylab.grid(True)
pylab.xlim(0.01,6000)
pylab.ylim(1e-9,1e-1)
pylab.ylabel(r'Rel. Int. Noise [1/Hz$^{1/2}$]')
pylab.legend(loc=1,prop={'size':10})

pylab.savefig('plots/OMC_RIN_noise.png',bbox_inches='tight')


fignum=fignum+1
pylab.figure(fignum)

pylab.subplot(2,1,1)

pylab.loglog(freq4,RIN_low,'r-',linewidth=0.8,label='OMC DCPD RIN, suppressed')
pylab.loglog(freq2,RIN,'b-',linewidth=0.8,label='OMC DCPD RIN, free-running')
pylab.loglog(freq8,dark,'-',color='0.6',linewidth=0.8,label=r'OMC DCPD Dark ($\div$ by SUM)')
pylab.loglog(freq2,NULL,'k-',linewidth=0.8,label=r'OMC DCPD NULL ($\div$ by SUM)')
pylab.loglog([min(freq2),max(freq2)],[5.6e-9, 5.6e-9],'y--',linewidth=1.4,label='Shot Noise, 10.2mA')
pylab.grid(True)
pylab.xlim(0.01,6000)
pylab.ylim(1e-9,1e-1)
pylab.ylabel(r'Rel. Int. Noise [1/Hz$^{1/2}$]')
pylab.legend(loc=1,prop={'size':8})

pylab.subplot(2,1,2)

pylab.semilogx(freq3,ISS,'b-',linewidth=0.8,label='ISS 2nd Loop PD5-8')
pylab.semilogx(freq3,POP,'g-',linewidth=0.8,label='POP A LF')
pylab.semilogx(freq3,ASAIR,'m-',linewidth=0.8,label='ASAIR A LF')
pylab.grid(True)
pylab.ylim(0,1)
pylab.xlim(0.01,6000)
pylab.ylabel('Coherence')
pylab.xlabel('Frequency [Hz]')
pylab.legend(loc=2,prop={'size':8})

pylab.savefig('plots/OMC_RIN_noise_coh.png',bbox_inches='tight')




fignum=fignum+1
pylab.figure(fignum)

pylab.loglog(freq1,deltaL,'b-',linewidth=0.8,label='OMC Length noise, free-running')
pylab.loglog(freq5,deltaL_low,'r-',linewidth=0.8,label='OMC Length noise, suppressed')
pylab.loglog([freq5[0],freq5[-1]],[3e-16,3e-16],'g--',linewidth=2)
pylab.grid(True)
pylab.xlim(0.01,7000)
pylab.ylim(8e-17,1e-9)
pylab.ylabel(r'Displacement noise [m / Hz$^{1/2}$]')
pylab.xlabel('Frequency [Hz]')
pylab.legend(loc=1,prop={'size':10})

pylab.savefig('plots/OMC_length_noise.png',bbox_inches='tight')




fignum=fignum+1
pylab.figure(fignum)

pylab.loglog(freq6,deltaL_L1,'r-',linewidth=0.8,label='L1 OMC, Sep 12 2014 (dither error signal, single bounce)')
pylab.loglog(freq7,deltaL_H1_PDH,'k-',linewidth=0.8,label='H1 OMC, Jan 13 2015 (PDH method, single bounce)')
pylab.loglog(freq1,deltaL,'b-',linewidth=0.8,label='H1 OMC, Jun 13 2015 (half-fringe method, full IFO)')

pylab.grid(True)
pylab.xlim(0.01,100000)
pylab.ylim(3e-17,1e-9)
pylab.ylabel(r'Displacement noise [m / Hz$^{1/2}$]')
pylab.xlabel('Frequency [Hz]')
pylab.legend(loc=1,prop={'size':10})

pylab.savefig('plots/OMC_length_noise_compare.png',bbox_inches='tight')


fignum=fignum+1
pylab.figure(fignum)

pylab.loglog(freq6,deltaL_L1,'r-',linewidth=0.8,label='L1 OMC, Sep 12 2014 (dither error signal, single bounce)')
pylab.loglog(freq7,deltaL_H1_PDH,'k-',linewidth=1.2,label='H1 OMC, Jan 13 2015 (PDH method, single bounce)')
pylab.loglog(freq1,deltaL,'b-',linewidth=1.2,label='H1 OMC, Jun 13 2015 (half-fringe method, full IFO)')

pylab.xscale('linear')
pylab.grid(True)
pylab.xlim(200,1800)
pylab.ylim(3e-17,1e-11)
pylab.ylabel(r'Displacement noise [m / Hz$^{1/2}$]')
pylab.xlabel('Frequency [Hz]')
pylab.legend(loc=1,prop={'size':10})

pylab.savefig('plots/OMC_length_noise_compare_zoom.png',bbox_inches='tight')



fignum=fignum+1
pylab.figure(fignum)

pylab.loglog(freq9,CAL,'k-',linewidth=0.8,label='Differential arm length sensitivity')
pylab.loglog(f,OMC_DARM,'r-',linewidth=0.8,label='OMC length noise contribution')

pylab.grid(True,which='major',linestyle=':')
pylab.grid(True,which='minor',linestyle=':')
pylab.xlim(5,6000)
pylab.ylim(1e-22,2e-15)
pylab.ylabel(r'Displacement noise [m / Hz$^{1/2}$]')
pylab.xlabel('Frequency [Hz]')
pylab.legend(loc=1,prop={'size':10})

pylab.savefig('plots/OMC_length_noise_project.png',bbox_inches='tight')



fignum=fignum+1
pylab.figure(fignum)

pylab.loglog(f,OMC_DARM,'r-',linewidth=0.8,label='OMC length noise --> DARM_ERR (Quadratic coupling)')
pylab.loglog(f_rms,OMC_DARM_rms,'b-',linewidth=0.8,label='OMC length noise --> DARM_ERR (Bilinear coupling)')
pylab.loglog(freq9,CAL,'k-',linewidth=0.8,label='H1:CAL-DELTAL_EXTERNAL')

pylab.grid(True,which='major',linestyle=':')
pylab.grid(True,which='minor',linestyle=':')
pylab.xlim(5,6000)
pylab.ylim(1e-22,2e-15)
pylab.ylabel(r'Displacement noise [m / Hz$^{1/2}$]')
pylab.xlabel('Frequency [Hz]')
pylab.legend(loc=1,prop={'size':10})

pylab.savefig('plots/OMC_length_noise_project_RMS.png',bbox_inches='tight')




OMC_NB = vstack((f,OMC_DARM)).T

savetxt('OMC_displacement_to_DARM_quadratic.txt',OMC_NB,fmt='%10.5e',delimiter='\t')


OMC_NB_rms = vstack((f_rms,OMC_DARM_rms)).T

savetxt('OMC_displacement_to_DARM_bilinear.txt',OMC_NB_rms,fmt='%10.5e',delimiter='\t')
