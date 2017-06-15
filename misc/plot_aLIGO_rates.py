#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab

def Nevents(range,rate,time):

    duty_cycle = 0.5
    return rate * (4*pi/3) * range**3 * time * duty_cycle / 1e6


IFO_sense = arange(10,300,1)

# Rates from LVC paper
# merger rates are in units of Mpc^-3 Myr^-1

lowrate = 0.01
midrate = 1.0
highrate = 10.0

# Fong et al rates: 90 - 270 - 1850 Gpc^-3 yr^-1
# https://arxiv.org/abs/1509.02922
# sGRB rates - seems likely that are counting extragalactic events? if some are SGR rate is too high
# Not all BNS might generate sGRB, so rate might be too low
#lowrate = (270-180) / 1e9 * 1e6
#midrate = 270 / 1e9 * 1e6
#highrate = (270+1580) / 1e9 * 1e6

fignum=0

fontP = FontProperties()
fontP.set_size('small')

fig_width_pt = 600  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio, = 0.618
golden_mean = 0.7         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width,fig_height]
matplotlib.rcParams.update(
        {'axes.labelsize': 14, 
            'font.size':   14, 
            'legend.fontsize': 14, 
            'xtick.labelsize': 14, 
            'ytick.labelsize': 14, 
            'text.usetex': False,
            'figure.figsize': fig_size,
            'font.family': "serif",
            'font.serif': ["Times New Roman"],
            'savefig.dpi': 250,
            'xtick.major.size':8,
            'xtick.minor.size':4,
            'ytick.major.size':8,
            'ytick.minor.size':4
            })

fignum+=1
pylab.figure(fignum)


# S5: 1yr triple, perhaps 0.5yr double? avg range ~15
# S6: 0.56yr 2+ detectors, avg range 18?
# call it 2.8 years to correct for assumption of 70% duty factor
# also say 17Mpc, some of this was triple coinc which would have better network SNR

r = 17
t = 2/0.7  # correction for duty factor assumption
y = Nevents(r,midrate,t)
yerr = array([[y-Nevents(r,lowrate,t)], [Nevents(r,highrate,t)-y]])
#yerr = [y-Nevents(r,lowrate,t), Nevents(r,highrate,t)-y]
pylab.errorbar(r, y, yerr=yerr, fmt='o--', ecolor='k', capthick=1)


"""
r = array([60.0, 70.0, 80.0])
t = 0.25
y = Nevents(r,midrate,t)
#yerr = array([[y-Nevents(r,lowrate,t)], [Nevents(r,highrate,t)-y]])
yerr = [y-Nevents(r,lowrate,t), Nevents(r,highrate,t)-y]
pylab.errorbar(r, y, yerr=yerr, fmt='o--', ecolor='k', capthick=1)
"""

r = 70
t = 0.27
y = Nevents(r,midrate,t)
yerr = array([[y-Nevents(r,lowrate,t)], [Nevents(r,highrate,t)-y]])
#yerr = [y-Nevents(r,lowrate,t), Nevents(r,highrate,t)-y]
pylab.errorbar(r, y, yerr=yerr, fmt='o--', ecolor='k', capthick=1)


r = array([80.0, 90.0, 100.0, 110.0, 120.0])
t = 0.5
y = Nevents(r,midrate,t)
#yerr = array([[y-Nevents(r,lowrate,t)], [Nevents(r,highrate,t)-y]])
yerr = [y-Nevents(r,lowrate,t), Nevents(r,highrate,t)-y]
pylab.errorbar(r, y, yerr=yerr, fmt='o--', ecolor='k', capthick=1)

r = array([130.0, 140.0, 150.0, 160.0, 170.0])
t = 0.75
y = Nevents(r,midrate,t)
#yerr = array([[y-Nevents(r,lowrate,t)], [Nevents(r,highrate,t)-y]])
yerr = [y-Nevents(r,lowrate,t), Nevents(r,highrate,t)-y]
pylab.errorbar(r, y, yerr=yerr, fmt='o--', ecolor='k', capthick=1)

r = 200
t = 1
y = Nevents(r,midrate,t)
yerr = array([[y-Nevents(r,lowrate,t)], [Nevents(r,highrate,t)-y]])
#yerr = [y-Nevents(r,lowrate,t), Nevents(r,highrate,t)-y]
pylab.errorbar(r, y, yerr=yerr, fmt='o--', ecolor='k', capthick=1)

pylab.grid(True, which='both', linestyle=':')

pylab.plot([10,220],[0.7,0.7],'k--',linewidth=1.1)

pylab.yscale('log')
pylab.ylim(3e-4,400)
pylab.xlim(10,220)

pylab.suptitle('Prospects for BNS detections, 50% duty cycle',fontsize=14)

#pylab.title('Error bars correspond to rates from Fong et al., arXiv:1509.02922. Except for S5+S6 a 2-detector network is assumed.\nRemember Poisson statistics, the probability of >= 1 detections when the expectation value = 0.7 is 50%.',fontsize=9,ha='center')

pylab.title('Note: error bars do NOT represent 67% CLs, they correspond to the "low" and "high" rates from arXiv:1003.2480.  Plot inspired by T1200307.\nExcept for S5+S6 a 2-detector network is assumed. Remember Poisson statistics, the probability of >= 1 detections when the expectation value = 0.7 is 50%.',fontsize=8,ha='center')

pylab.xlabel('Average IFO Range [Mpc]')
pylab.ylabel('Expected Number of Events')

pylab.legend(('S5+S6','O1 (4 months)','6 months','9 months','1 year'),loc=2,prop=fontP,fancybox=True)

pylab.savefig('Nevents_aLIGO.pdf')
pylab.savefig('Nevents_aLIGO.png')
