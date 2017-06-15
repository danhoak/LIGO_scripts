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

data1 = genfromtxt('data/OMC_OLG_offset_13Jun2015_A.txt')
freq1 = data1[4:,0]
mag1 = data1[4:,1]
phase1 = data1[4:,2]*pi/180


data2 = genfromtxt('data/OMC_OLG_offset_13Jun2015_B.txt')
freq2 = data2[:,0]
mag2 = data2[:,1]
phase2 = data2[:,2]*pi/180
 
freq = hstack((freq2,freq1))
mag = hstack((mag2,mag1))
phase = hstack((phase2,phase1))

data = vstack((data2[:,:],data1[4:,:]))


fignum=0

fignum=fignum+1
pylab.figure(fignum)

pylab.subplot(2,1,1)

pylab.loglog(freq,mag,'bo-')
pylab.grid(True)
pylab.xlim(0.6,2000)
pylab.ylabel('Magnitude')


pylab.subplot(2,1,2)

pylab.semilogx(freq,phase*180/pi,'bo-')
pylab.grid(True)
pylab.ylim(-180,180)
pylab.xlim(0.6,2000)
pylab.ylabel('Phase [deg]')
pylab.xlabel('Frequency [Hz]')

pylab.savefig('figures/OMC_OLTF.png')

print shape(data)

CLG_inv = ones(shape(data))
CLG_inv_meters = ones(shape(data))

RINfactor = 2.72e-9 # meters per RIN

for i in range(len(freq)):

    # calculate the denominator of the CLG
    CLG = 1 + mag[i]*exp(1j*phase[i])
    CLG_mag = abs(CLG)
    CLG_phase = cm.phase(CLG)

    #print freq[i], CLG_mag, CLG_phase*180/pi

    CLG_inv[i,0] = freq[i]
    CLG_inv_meters[i,0] = freq[i]

    #CLG_inv_meters[i,1] = 20*log10(CLG_mag*RINfactor) # mag in dB, conv to meters
    CLG_inv_meters[i,1] = CLG_mag*RINfactor # mag in dB, conv to meters

    #CLG_inv[i,1] = 20*log10(CLG_mag) # mag in dB, conv to RIN
    CLG_inv[i,1] = CLG_mag # mag in...uh, mag

    # We should only need the magnitude information, but carry the phase along for completeness

    CLG_inv[i,2] = CLG_phase
    CLG_inv_meters[i,2] = CLG_phase

savetxt('data/OMC_CLG_13Jun_RIN.txt', CLG_inv, fmt='%10.5f', delimiter='   ')
savetxt('data/OMC_CLG_13Jun_meters.txt', CLG_inv_meters, fmt='%10.5e', delimiter='   ')

