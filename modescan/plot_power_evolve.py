#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
#from scipy import signal
#import matplotlib.gridspec as gridspec
matplotlib.rcParams.update({'savefig.dpi':250})

x = genfromtxt('data/low2hi_data.txt')
t = x[:,0]
SB9_power = x[:,1]
SB45_power = x[:,2]
CR_power = x[:,3]


fignum=0
fignum=fignum+1
pylab.figure(fignum)

pylab.subplot(3,1,1)
pylab.plot(t,CR_power,'bo',markersize=5)
pylab.plot([-50,-50],[-1,18],'r--',linewidth=1.4)
pylab.plot([-150,-150],[-1,18],'r--',linewidth=1.4)
pylab.ylim(1.4,2.0)
pylab.grid(True, which='both', linestyle=':', alpha=0.4)
pylab.ylabel('Carrier', fontsize=10)
pylab.xticks(visible=False)
pylab.yticks(fontsize=8)
pylab.title('Power content at OMC, 2.2W vs 22.6W (t=0 is end of power-up)\nData are DCPD photocurrent normalized to input power',fontsize=12)

pylab.subplot(3,1,2)
pylab.plot(t,SB45_power,'bo',markersize=5)
pylab.plot([-50,-50],[6,18],'r--',linewidth=1.4)
pylab.plot([-150,-150],[6,18],'r--',linewidth=1.4)
pylab.ylim(6,18)
pylab.arrow(-180, 12, -350, 0, head_width=0.5, head_length=50, fc='r', ec='r')
pylab.arrow(-20, 12, 420, 0, head_width=0.5, head_length=50, fc='r', ec='r')
pylab.annotate('2.2W', xy = (-230,13.5), color='r', xytext = (1, 0), textcoords = 'offset points', ha = 'right', va = 'top',fontsize=10)
pylab.annotate('22.6W', xy = (340,13.5), color='r', xytext = (1, 0), textcoords = 'offset points', ha = 'right', va = 'top',fontsize=10)
pylab.ylabel('45MHz Sidebands', fontsize=10)
pylab.grid(True, which='both', linestyle=':', alpha=0.4)
pylab.xticks(visible=False)
pylab.yticks(fontsize=8)

pylab.subplot(3,1,3)
pylab.plot(t,SB9_power,'bo',markersize=5)
pylab.plot([-50,-50],[0,18],'r--',linewidth=1.4)
pylab.plot([-150,-150],[0,18],'r--',linewidth=1.4)
pylab.ylim(0.2,1.7)
pylab.grid(True, which='both', linestyle=':', alpha=0.4)
pylab.ylabel('9MHz Sidebands', fontsize=10)
pylab.xticks(fontsize=8)
pylab.yticks(fontsize=8)
pylab.xlabel('Sweep time [sec]',fontsize=10)

pylab.savefig('plots/low2hi_power_data.png')
