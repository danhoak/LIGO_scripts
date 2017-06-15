#! /usr/bin/env python
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
from numpy import *
import matplotlib
from matplotlib.font_manager import FontProperties
matplotlib.use("Agg")
import pylab
from scipy.special import jn

fig_width_pt = 600  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
#golden_mean = (2.236-1.0)/2.0         # Aesthetic ratio, = 0.618
golden_mean = 0.6
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size = [fig_width,fig_height]
matplotlib.rcParams.update({'figure.figsize': fig_size})

c = 299792458.0

#dx = 0.01e-10
#x = arange(-5e-9,5e-9,dx)
#lam = 1064.0e-9

dx = 0.0001
x = arange(-0.5,0.5,dx)
lam = 1064.0/10

w = 0.02
s = -0.079
x1 = arange(s,s+w,dx)
x2 = arange(-w/2,w/2,dx)

E_in = 1.0

# OMC input & output couplers are 8000ppm transmission

r_i = sqrt(1-8000e-6)
r_e = sqrt(1-8000e-6)

t_i = sqrt(1-r_i**2)
t_e = sqrt(1-r_e**2)


L = 0

finesse = pi*sqrt(r_i * r_e) / (1- r_i*r_e)

trans_pow = (t_i * t_e)**2 / (1- r_i*r_e)**2
E_trans = E_in * ( t_i * t_e * exp(1.0j*2*pi*(L+x)/lam))/(1 - r_i * r_e * exp(2.0j*2*pi*(L+x)/lam))
P_trans = real(E_trans * conj(E_trans))/trans_pow

E_trans1 = E_in * ( t_i * t_e * exp(1.0j*2*pi*(L+x1)/lam))/(1 - r_i * r_e * exp(2.0j*2*pi*(L+x1)/lam))
P_trans1 = real(E_trans1 * conj(E_trans1))/trans_pow

E_trans2 = E_in * ( t_i * t_e * exp(1.0j*2*pi*(L+x2)/lam))/(1 - r_i * r_e * exp(2.0j*2*pi*(L+x2)/lam))
P_trans2 = real(E_trans2 * conj(E_trans2))/trans_pow

dPdx = gradient(P_trans)

print
print 'Cavity finesse is %.3f' % finesse

halfpoint = arccos(1 - (1-r_i*r_e)**2 / (2*r_i*r_e))
x_2 = lam*halfpoint/(4*pi)

# Calculate the halfway point in two ways: first invert the transmission expression, second use the finesse (?)

print('halfway off-resonance [m] = %g' % x_2)

xr = lam/(4*finesse)

print('alternate for halfway [m] = %g' % xr)

idx = (abs(x-x_2)).argmin()
ridx = (abs(x-xr)).argmin()

print('halfway off-resonance trans: %g' % (P_trans[idx]))
print('halfway off-resonance optical gain [RIN/m]: %g' % (-0.5*dPdx[idx]/dx))
#print('halfway off-resonance optical gain [RIN/m]: %g' % (-0.5*dPdx[ridx]/dx))

print('meters/RIN for offset lock: %g' % (1/((-0.5*dPdx[idx]/dx))))
print

# measuring RIN -- need to correct for 50% transmission at halfway point?
# check whitening gain, transimpedance gain...measure whitening gain?

matplotlib.rcParams.update({'savefig.dpi':250})

fignum=0

fignum=fignum+1
pylab.figure(fignum)

pylab.plot(x,P_trans,'k-')
pylab.plot(x1,P_trans1,'k-',linewidth=3)
pylab.plot(x2,P_trans2,'k-',linewidth=3)
#pylab.plot(x,P_trans2,'r-',linewidth=2.0)
#pylab.plot([x2,x2],[0,1],'r--')

dx = x1[0]-x1[1]
dy = P_trans1[0]-P_trans1[1]
pylab.arrow(x1[1], P_trans1[1], dx, dy, head_width=0.015, head_length=0.025, fc='k', ec='k')

dx = x1[-2]-x1[-1]
dy = P_trans1[-2]-P_trans1[-1]
pylab.arrow(x1[-2], P_trans1[-2], -dx, -dy, head_width=0.015, head_length=0.025, fc='k', ec='k')


dx = x2[0]-x2[1]
dy = P_trans2[0]-P_trans2[1]
pylab.arrow(x2[1], P_trans2[1], dx, dy, head_width=0.015, head_length=0.025, fc='k', ec='k')

dx = x2[-2]-x2[-1]
dy = P_trans2[-2]-P_trans2[-1]
pylab.arrow(x2[-2], P_trans2[-2], -dx, -dy, head_width=0.015, head_length=0.025, fc='k', ec='k')


pylab.plot([x1[0], x1[0]],[0.0, P_trans1[0]],'k:')
pylab.plot([x1[-1], x1[-1]],[0.0, P_trans1[-1]],'k:')

pylab.plot([x2[0], x2[0]],[0.0, P_trans2[0]],'k:')
pylab.plot([x2[-1], x2[-1]],[0.0, P_trans2[-1]],'k:')

pylab.plot([x[0], x1[0]],[P_trans1[0], P_trans1[0]],'k:')
pylab.plot([x[0], x1[-1]],[P_trans1[-1], P_trans1[-1]],'k:')

pylab.plot([x[0], x2[0]],[min(P_trans2), min(P_trans2)],'k:')
pylab.plot([x[0], x2[-1]],[max(P_trans2), max(P_trans2)],'k:')

#pylab.grid(True,which='both', linestyle='--')
pylab.xlim(-0.4,0.4)
#pylab.yscale('log')
pylab.ylim(0.0,1.05)

pylab.text(x1[0]-0.0025,-0.05,r'$\delta L$',fontsize=16)
pylab.text(x2[0]-0.0025,-0.05,r'$\delta L$',fontsize=16)

pylab.text(-0.45,mean(P_trans1)-0.01,r'$\delta P$',fontsize=16)
pylab.text(-0.45,0.975,r'$\delta P$',fontsize=16)

pylab.xticks(visible=False)
pylab.yticks(visible=False)
#pylab.xtick_labels(visible=False)

pylab.ticklabel_format(style = 'sci', useOffset=False)

pylab.savefig('plots/cavity_offset_RIN.png',bbox_inches='tight')
